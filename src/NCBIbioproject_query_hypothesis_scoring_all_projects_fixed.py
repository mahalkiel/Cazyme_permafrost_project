#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 09:32:00 2026
Purpose: Query NCBI BioProject for candidate Arctic/subarctic metagenomic
datasets and retrieve project-level metadata.

@author: michaelsnyder, Codex assisted
"""

from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
import time
import re
import os
from collections import defaultdict
from pathlib import Path

from eco_harmonization_integrated import apply_harmonization

# ---------------------------------------------------------------------------
# NCBI / Entrez settings
# ---------------------------------------------------------------------------

Entrez.email = "michaelsnyder777@gmail.com"

# Optional:
Entrez.api_key = os.getenv("NCBI_API_KEY")

DB_NAME = "bioproject"
SEARCH_PAGE_SIZE = 500  # Batch size only; no overall BioProject limit
REQUEST_DELAY = 0.34
BATCH_SIZE = 200
EXCLUDED_LIBRARY_STRATEGIES = {"AMPLICON"}

OUTPUT_FOLDER = Path(__file__).resolve().parent.parent / "output" / "explr_metadata"
print(OUTPUT_FOLDER)

OUTPUT_PROJECTS = OUTPUT_FOLDER / "candidate_bioprojects.csv"
OUTPUT_BIOSAMPLES = OUTPUT_FOLDER / "candidate_biosamples.csv"
OUTPUT_ATTRIBUTES = OUTPUT_FOLDER / "candidate_biosample_attributes.csv"
OUTPUT_SRA = OUTPUT_FOLDER / "candidate_sra_runs.csv"
OUTPUT_SCREENING = OUTPUT_FOLDER / "candidate_cazyme_screening.csv"
OUTPUT_HARMONIZED = OUTPUT_FOLDER / "candidate_cazyme_screening_harmonized.csv"
OUTPUT_CANDIDATE_REVIEW = OUTPUT_FOLDER / "candidate_review.csv"


# ---------------------------------------------------------------------------
# Exploratory BioProject query
# ---------------------------------------------------------------------------

SEARCH_TERM = """\
(
  permafrost[All Fields]
  OR arctic[All Fields]
  OR subarctic[All Fields]
  OR tundra[All Fields]
  OR thermokarst[All Fields]
  OR cryosol[All Fields]
)
AND
(
  metagenome[All Fields]
  OR metagenomic[All Fields]
)
NOT
(
  gut[All Fields]
  OR fecal[All Fields]
  OR human[All Fields]
  OR marine[All Fields]
  OR seawater[All Fields]
  OR sea water[All Fields]
  OR oceanic[All Fields]
  OR meltwater[All Fields]
  OR melt water[All Fields]
  OR sea ice[All Fields]
  OR sea-ice[All Fields]
  OR melted sea ice[All Fields]
  OR melted sea-ice[All Fields]
  OR metatranscriptome[All Fields]
  OR metatranscriptomic[All Fields]
  OR transcriptome[All Fields]
  OR transcriptomic[All Fields]
  OR snow[All Fields]
  OR Pacific Ocean[All Fields]
)"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def text_at(element, xpath):
    node = element.find(xpath)
    if node is not None and node.text:
        return node.text.strip()
    return None


def combined_text(element):
    if element is None:
        return None

    value = " ".join(
        x.strip()
        for x in element.itertext()
        if x and x.strip()
    )
    return value or None


def search_bioprojects(search_term, page_size=SEARCH_PAGE_SIZE):
    """Return every matching BioProject UID, paginating the Entrez history set.

    page_size controls the size of each request, NOT the total projects screened.
    No SRA sequence reads are downloaded by this function or pipeline.
    """
    print(f"Searching NCBI {DB_NAME} database...")
    with Entrez.esearch(db=DB_NAME, term=search_term,
                        retmax=0, usehistory="y") as handle:
        record = Entrez.read(handle)

    total = int(record.get("Count", 0))
    print(f"Total query matches: {total}")
    if not total:
        return []

    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    uids = []
    for start in range(0, total, page_size):
        count = min(page_size, total - start)
        with Entrez.esearch(db=DB_NAME, term=f"#{query_key}",
                            WebEnv=webenv, retstart=start, retmax=count,
                            usehistory="y") as handle:
            page = Entrez.read(handle)
        batch = list(page.get("IdList", []))
        if len(batch) != count:
            raise RuntimeError(
                f"Incomplete Entrez page at offset {start}: "
                f"expected {count} UIDs, got {len(batch)}"
            )
        uids.extend(batch)
        print(f"Retrieved {len(uids)}/{total} BioProject UIDs")
        time.sleep(REQUEST_DELAY)

    if len(set(uids)) != total:
        raise RuntimeError("BioProject UID pagination produced duplicates or missing records")
    return uids


def fetch_bioproject_xml(uid):
    """Fetch a complete BioProject XML record by internal NCBI UID."""
    with Entrez.efetch(
        db=DB_NAME,
        id=str(uid),
        rettype="xml",
        retmode="xml"
    ) as handle:
        data = handle.read()

    return ET.fromstring(data)


def parse_bioproject_xml(root, uid):
    """
    Parse project-level BioProject metadata.

    Fields:
      - accession
      - internal NCBI UID
      - title
      - description
      - organism
      - environment
      - project type
      - submission date
      - linked umbrella / parent projects
    """

    project = root.find(".//Project")

    if project is None:
        raise ValueError(f"No <Project> element found for UID {uid}")

    # Public BioProject accession
    archive = project.find("./ProjectID/ArchiveID")

    accession = None
    archive_id = None

    if archive is not None:
        accession = archive.attrib.get("accession")
        archive_id = archive.attrib.get("id")

    title = text_at(project, "./ProjectDescr/Title")
    description = text_at(project, "./ProjectDescr/Description")

    # Organism
    organism = (
        text_at(project, ".//Organism/OrganismName")
        or text_at(project, ".//OrganismName")
    )

    taxonomy_id = None
    organism_node = project.find(".//Organism")

    if organism_node is not None:
        taxonomy_id = (
            organism_node.attrib.get("taxID")
            or organism_node.attrib.get("taxonomy_id")
        )

    # Environment
    environment = None

    for xpath in (
        ".//EnvironmentalSample",
        ".//Environment",
        ".//Target",
    ):
        node = project.find(xpath)
        value = combined_text(node)

        if value:
            environment = value
            break

    # Project type
    project_type_node = project.find("./ProjectType")
    project_type_parts = []

    if project_type_node is not None:
        for child in project_type_node.iter():

            if child is project_type_node:
                continue

            tag = child.tag.split("}")[-1]
            text = (child.text or "").strip()

            if len(child) == 0:
                value = f"{tag}: {text}" if text else tag

                if value not in project_type_parts:
                    project_type_parts.append(value)

    project_type = (
        " | ".join(project_type_parts)
        if project_type_parts
        else None
    )

    # Submission date
    submission_date = None

    for node in root.iter():
        tag = node.tag.split("}")[-1].lower()

        if tag == "submission":
            for attr in (
                "submitted",
                "submission_date",
                "date"
            ):
                if node.attrib.get(attr):
                    submission_date = node.attrib[attr]
                    break

        if submission_date:
            break

    if submission_date is None:
        for xpath in (
            ".//SubmissionDate",
            ".//Submission/SubmissionDate",
        ):
            submission_date = text_at(root, xpath)

            if submission_date:
                break

    # Umbrella / parent projects
    umbrella_projects = []

    for parent_tag in (
        "ParentProject",
        "UmbrellaProject",
        "SuperProject",
    ):
        for parent in project.findall(f".//{parent_tag}"):

            for descendant in parent.iter():
                candidate = descendant.attrib.get("accession")

                if (
                    candidate
                    and candidate != accession
                    and candidate not in umbrella_projects
                ):
                    umbrella_projects.append(candidate)

    # Additional fallback for parent ArchiveID nodes
    for node in project.findall(".//ArchiveID"):
        candidate = node.attrib.get("accession")

        attributes = " ".join(
            str(value).lower()
            for value in node.attrib.values()
        )

        if (
            candidate
            and candidate != accession
            and any(
                token in attributes
                for token in ("parent", "umbrella", "super")
            )
            and candidate not in umbrella_projects
        ):
            umbrella_projects.append(candidate)

    return {
        "bioproject_accession": accession,
        "bioproject_uid": str(uid),
        "archive_id": archive_id,
        "title": title,
        "description": description,
        "organism": organism,
        "taxonomy_id": taxonomy_id,
        "environment": environment,
        "project_type": project_type,
        "submission_date": submission_date,
        "umbrella_projects": "; ".join(umbrella_projects) or None,
    }


def retrieve_candidate_metadata(uids):
    """Fetch and parse project metadata for all candidate BioProject UIDs."""

    rows = []

    for index, uid in enumerate(uids, start=1):

        print(
            f"[{index}/{len(uids)}] "
            f"Fetching BioProject UID {uid}"
        )

        try:
            root = fetch_bioproject_xml(uid)

            row = parse_bioproject_xml(
                root=root,
                uid=uid
            )

            row["retrieval_error"] = None
            rows.append(row)

        except Exception as exc:

            print(
                f"  Warning: could not process "
                f"UID {uid}: {exc}"
            )

            rows.append({
                "bioproject_accession": None,
                "bioproject_uid": str(uid),
                "archive_id": None,
                "title": None,
                "description": None,
                "organism": None,
                "taxonomy_id": None,
                "environment": None,
                "project_type": None,
                "submission_date": None,
                "umbrella_projects": None,
                "retrieval_error": str(exc),
            })

        time.sleep(REQUEST_DELAY)

    return pd.DataFrame(rows)



# ---------------------------------------------------------------------------
# BioProject -> BioSample / SRA retrieval
# ---------------------------------------------------------------------------

def chunked(values, size=BATCH_SIZE):
    """Yield successive fixed-size batches from a list."""
    for start in range(0, len(values), size):
        yield values[start:start + size]


def linked_ids(dbfrom, db, source_id):
    """Return NCBI IDs linked from one Entrez database to another."""
    with Entrez.elink(dbfrom=dbfrom, db=db, id=str(source_id)) as handle:
        result = Entrez.read(handle)

    ids = []
    for linkset in result:
        for link_db in linkset.get("LinkSetDb", []):
            for link in link_db.get("Link", []):
                linked_id = link.get("Id")
                if linked_id:
                    ids.append(str(linked_id))

    return list(dict.fromkeys(ids))


def fetch_xml_batch(db, ids):
    """Fetch an XML batch from an Entrez database."""
    if not ids:
        return None

    with Entrez.efetch(
        db=db,
        id=",".join(map(str, ids)),
        rettype="xml",
        retmode="xml"
    ) as handle:
        data = handle.read()

    return ET.fromstring(data)


def clean_column_name(value):
    """Create a stable dataframe-friendly name for BioSample attributes."""
    value = str(value).strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_") or "unnamed_attribute"


def parse_biosample_batch(root, bioproject_accession, bioproject_uid):
    """
    Parse BioSample XML without discarding submitted environmental attributes.

    Returns
    -------
    sample_rows : list[dict]
        One row per BioSample, with all attributes widened into columns.
    attribute_rows : list[dict]
        Lossless long-form representation of every submitted attribute.
    """
    sample_rows = []
    attribute_rows = []

    for sample in root.findall(".//BioSample"):
        accession = sample.attrib.get("accession")
        sample_id = sample.attrib.get("id")

        row = {
            "bioproject_accession": bioproject_accession,
            "bioproject_uid": str(bioproject_uid),
            "biosample_accession": accession,
            "biosample_uid": sample_id,
            "sample_title": text_at(sample, "./Description/Title"),
            "sample_comment": combined_text(sample.find("./Description/Comment")),
            "organism": text_at(sample, "./Description/Organism/OrganismName"),
            "taxonomy_id": None,
            "package": None,
        }

        organism_node = sample.find("./Description/Organism")
        if organism_node is not None:
            row["taxonomy_id"] = (
                organism_node.attrib.get("taxonomy_id")
                or organism_node.attrib.get("taxid")
                or organism_node.attrib.get("taxID")
            )

        models = sample.find("./Models")
        if models is not None:
            model_names = [
                (model.text or "").strip()
                for model in models.findall("./Model")
                if (model.text or "").strip()
            ]
            row["package"] = "; ".join(model_names) or None

        # Preserve every submitted/harmonized BioSample attribute.
        seen_names = defaultdict(int)
        for attribute in sample.findall(".//Attributes/Attribute"):
            submitted_name = attribute.attrib.get("attribute_name")
            harmonized_name = attribute.attrib.get("harmonized_name")
            display_name = harmonized_name or submitted_name or "unnamed_attribute"
            value = (attribute.text or "").strip() or None
            unit = attribute.attrib.get("unit")

            base_column = clean_column_name(display_name)
            seen_names[base_column] += 1
            column = base_column
            if seen_names[base_column] > 1:
                column = f"{base_column}_{seen_names[base_column]}"

            row[column] = value

            attribute_rows.append({
                "bioproject_accession": bioproject_accession,
                "bioproject_uid": str(bioproject_uid),
                "biosample_accession": accession,
                "biosample_uid": sample_id,
                "attribute_name": submitted_name,
                "harmonized_name": harmonized_name,
                "normalized_name": base_column,
                "value": value,
                "unit": unit,
            })

        sample_rows.append(row)

    return sample_rows, attribute_rows


def retrieve_biosamples(projects_df):
    """Retrieve and parse linked BioSamples for each successfully parsed project."""
    sample_rows = []
    attribute_rows = []

    valid = projects_df.dropna(subset=["bioproject_accession", "bioproject_uid"])

    for index, project in enumerate(valid.itertuples(index=False), start=1):
        accession = project.bioproject_accession
        uid = project.bioproject_uid
        print(f"[{index}/{len(valid)}] Linking {accession} -> BioSample")

        try:
            ids = linked_ids("bioproject", "biosample", uid)
            print(f"  Linked BioSamples: {len(ids)}")

            for batch in chunked(ids):
                root = fetch_xml_batch("biosample", batch)
                if root is None:
                    continue
                rows, attrs = parse_biosample_batch(root, accession, uid)
                sample_rows.extend(rows)
                attribute_rows.extend(attrs)
                time.sleep(REQUEST_DELAY)

        except Exception as exc:
            print(f"  Warning: BioSample retrieval failed for {accession}: {exc}")

        time.sleep(REQUEST_DELAY)

    samples_df = pd.DataFrame(sample_rows)
    attributes_df = pd.DataFrame(attribute_rows)
    return samples_df, attributes_df


def first_text(element, xpaths):
    """Return the first populated text value among several XPath candidates."""
    for xpath in xpaths:
        value = text_at(element, xpath)
        if value:
            return value
    return None


def parse_sra_batch(root, bioproject_accession, bioproject_uid):
    """Parse SRA EXPERIMENT_PACKAGE XML to one row per sequencing run."""
    rows = []

    for package in root.findall(".//EXPERIMENT_PACKAGE"):
        experiment = package.find("./EXPERIMENT")
        study_ref = package.find(".//STUDY_REF")
        sample_ref = package.find(".//SAMPLE_REF")

        common = {
            "bioproject_accession": bioproject_accession,
            "bioproject_uid": str(bioproject_uid),
            "study_accession": study_ref.attrib.get("accession") if study_ref is not None else None,
            "experiment_accession": experiment.attrib.get("accession") if experiment is not None else None,
            "biosample_accession": None,
            "sra_sample_accession": sample_ref.attrib.get("accession") if sample_ref is not None else None,
            "experiment_title": first_text(package, ["./EXPERIMENT/TITLE", ".//EXPERIMENT/TITLE"]),
            "library_name": first_text(package, [".//LIBRARY_DESCRIPTOR/LIBRARY_NAME"]),
            "library_strategy": first_text(package, [".//LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY"]),
            "library_source": first_text(package, [".//LIBRARY_DESCRIPTOR/LIBRARY_SOURCE"]),
            "library_selection": first_text(package, [".//LIBRARY_DESCRIPTOR/LIBRARY_SELECTION"]),
            "library_layout": None,
            "platform": None,
            "instrument_model": None,
        }

        # BioSample accession is often represented as an external ID.
        for external_id in package.findall(".//SAMPLE/IDENTIFIERS/EXTERNAL_ID"):
            namespace = (external_id.attrib.get("namespace") or "").lower()
            value = (external_id.text or "").strip()
            if value and ("biosample" in namespace or value.startswith("SAM")):
                common["biosample_accession"] = value
                break

        if not common["biosample_accession"]:
            for primary_id in package.findall(".//SAMPLE/IDENTIFIERS/PRIMARY_ID"):
                value = (primary_id.text or "").strip()
                if value.startswith("SAM"):
                    common["biosample_accession"] = value
                    break

        layout = package.find(".//LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT")
        if layout is not None and len(layout):
            common["library_layout"] = layout[0].tag.split("}")[-1]

        platform = package.find(".//PLATFORM")
        if platform is not None and len(platform):
            common["platform"] = platform[0].tag.split("}")[-1]
            common["instrument_model"] = text_at(platform, ".//INSTRUMENT_MODEL")

        runs = package.findall(".//RUN_SET/RUN")
        if not runs:
            runs = package.findall(".//RUN")

        for run in runs:
            row = common.copy()
            row.update({
                "run_accession": run.attrib.get("accession"),
                "total_spots": run.attrib.get("total_spots"),
                "total_bases": run.attrib.get("total_bases"),
                "run_size": run.attrib.get("size"),
                "published": run.attrib.get("published"),
            })
            rows.append(row)

    return rows


def retrieve_sra(projects_df):
    """Retrieve SRA records linked to each candidate BioProject."""
    rows = []
    valid = projects_df.dropna(subset=["bioproject_accession", "bioproject_uid"])

    for index, project in enumerate(valid.itertuples(index=False), start=1):
        accession = project.bioproject_accession
        uid = project.bioproject_uid
        print(f"[{index}/{len(valid)}] Linking {accession} -> SRA")

        try:
            ids = linked_ids("bioproject", "sra", uid)
            print(f"  Linked SRA records: {len(ids)}")

            for batch in chunked(ids):
                root = fetch_xml_batch("sra", batch)
                if root is None:
                    continue
                rows.extend(parse_sra_batch(root, accession, uid))
                time.sleep(REQUEST_DELAY)

        except Exception as exc:
            print(f"  Warning: SRA retrieval failed for {accession}: {exc}")

        time.sleep(REQUEST_DELAY)

    return pd.DataFrame(rows)


def filter_sra_records(sra_df):
    """Remove linked SRA records that are outside the CAZyme screening design."""
    if sra_df.empty or "library_strategy" not in sra_df.columns:
        return sra_df

    strategy = sra_df["library_strategy"].fillna("").str.upper()
    excluded = strategy.isin(EXCLUDED_LIBRARY_STRATEGIES)

    if excluded.any():
        counts = strategy[excluded].value_counts().to_dict()
        print(f"Excluded {int(excluded.sum())} SRA run(s) by library_strategy: {counts}")

    return sra_df.loc[~excluded].copy()


def build_screening_table(projects_df, biosamples_df, sra_df):
    """
    Create one screening row per SRA run, carrying project and sample metadata.

    This is intentionally a discovery table rather than a final harmonized
    ecological dataset. Raw BioSample attribute names are retained.
    """
    if sra_df.empty:
        return pd.DataFrame()

    project_cols = [
        "bioproject_accession", "title", "description", "organism",
        "environment", "project_type", "submission_date", "umbrella_projects"
    ]
    project_lookup = projects_df[project_cols].rename(columns={
        "title": "project_title",
        "description": "project_description",
        "organism": "project_organism",
    })

    screening = sra_df.merge(project_lookup, on="bioproject_accession", how="left")

    if not biosamples_df.empty and "biosample_accession" in biosamples_df.columns:
        sample_lookup = biosamples_df.drop_duplicates(
            subset=["bioproject_accession", "biosample_accession"]
        )
        screening = screening.merge(
            sample_lookup,
            on=["bioproject_accession", "biosample_accession"],
            how="left",
            suffixes=("", "_biosample")
        )

    # Lightweight technical screening flags; ecological harmonization comes later.
    strategy = screening.get("library_strategy", pd.Series(index=screening.index, dtype="object"))
    source = screening.get("library_source", pd.Series(index=screening.index, dtype="object"))
    screening["is_wgs"] = strategy.fillna("").str.upper().eq("WGS")
    screening["is_metagenomic"] = source.fillna("").str.upper().eq("METAGENOMIC")
    screening["cazyme_technical_candidate"] = screening["is_wgs"] & screening["is_metagenomic"]

    return screening



# ---------------------------------------------------------------------------
# Candidate-review schema and conservative biological-state screening
# ---------------------------------------------------------------------------

CANDIDATE_REVIEW_COLUMNS = [
    "bioproject_accession", "project_title", "biosample_accession", "sample_name",
    "run_accession", "experiment_accession", "library_strategy", "library_source",
    "platform", "instrument_model", "layout", "read_count", "base_count",
    "site_name", "country_region", "latitude", "longitude", "collection_date",
    "depth_raw", "depth_cm_min", "depth_cm_max", "ecosystem", "soil_type",
    "state_raw", "permafrost_state", "analysis_state", "state_source",
    "state_evidence", "state_confidence", "temperature_c", "ph",
    "moisture_percent", "site_comparison_states", "comparison_coverage",
    "depth_comparability", "candidate_score", "eligibility_tier",
    "exclusion_reason", "manual_review", "review_notes",
]

def first_nonempty(row, names):
    for name in names:
        if name not in row.index:
            continue
        value = row.get(name)
        if pd.notna(value) and str(value).strip() not in {"", "nan", "None"}:
            return value, name
    return None, None

def numeric_value(value):
    if value is None or pd.isna(value):
        return None
    match = re.search(r"-?\d+(?:\.\d+)?", str(value).replace(",", ""))
    return float(match.group()) if match else None

def parse_depth_for_review(value):
    if value is None or pd.isna(value):
        return None, None
    text = str(value).strip().lower().replace("_", "-")
    if not text:
        return None, None
    if text == "surface":
        return 0.0, 0.0
    m = re.search(r"(?P<a>\d+(?:\.\d+)?)\s*(?:-|to)\s*(?P<b>\d+(?:\.\d+)?)\s*(?P<u>cm|m|mm)?", text)
    if m:
        a, b = float(m.group("a")), float(m.group("b"))
        factor = {"mm": 0.1, "cm": 1.0, "m": 100.0}[m.group("u") or "cm"]
        return a * factor, b * factor
    m = re.search(r"(?P<x>\d+(?:\.\d+)?)\s*(?P<u>cm|m|mm)\b", text)
    if m:
        x = float(m.group("x")); factor = {"mm": 0.1, "cm": 1.0, "m": 100.0}[m.group("u")]
        return x * factor, x * factor
    return None, None

def classify_permafrost_state(row):
    """Infer permafrost state from sample-level evidence before project context.

    Project-level wording is deliberately contextual only: a project mentioning
    thermokarst, active layer, or permafrost may contain several ecological
    states and must not blanket-label every linked BioSample.
    """
    sample_fields = [
        "permafrost_state", "sample_type", "env_local_scale", "env_medium",
        "isolation_source", "sample_title", "sample_comment", "experiment_title",
        "treatment", "samp_mat_process",
    ]
    project_fields = ["project_title", "project_description"]

    sample_pieces = []
    for field in sample_fields:
        if field in row.index:
            value = row.get(field)
            if pd.notna(value) and str(value).strip():
                sample_pieces.append((field, str(value).strip()))
    sample_text = " | ".join(v.lower() for _, v in sample_pieces)

    project_pieces = []
    for field in project_fields:
        if field in row.index:
            value = row.get(field)
            if pd.notna(value) and str(value).strip():
                project_pieces.append((field, str(value).strip()))
    project_text = " | ".join(v.lower() for _, v in project_pieces)

    # Ordered from specific sample-state language to broader terminology.
    direct_rules = [
        ("transition_zone", ["transition zone", "transition layer", "transient layer"]),
        ("active_layer", ["active layer", "active-layer", "seasonally thawed layer"]),
        ("recently_thawed", ["recently thawed", "newly thawed", "abrupt thaw"]),
        ("degrading_permafrost", ["degrading permafrost", "thawing permafrost", "degraded permafrost"]),
        ("thermokarst", ["thermokarst"]),
        ("intact_permafrost", ["intact permafrost", "undisturbed permafrost", "frozen permafrost", "permafrost"]),
    ]
    for state, terms in direct_rules:
        matched = next((term for term in terms if term in sample_text), None)
        if matched:
            field, value = next(((f, v) for f, v in sample_pieces if matched in v.lower()), (None, None))
            analysis_state = (
                "intact_permafrost" if state == "intact_permafrost" else
                "active_layer" if state == "active_layer" else
                "thaw_affected" if state in {"degrading_permafrost", "recently_thawed", "thermokarst"} else
                "exclude"
            )
            return state, analysis_state, "NCBI_BioSample", f"{field}={value}", "direct"

    # Use harmonized soil-layer calls only after raw sample-state text fails.
    # These are inferred labels derived from BioSample metadata, so confidence
    # stays below direct submitter state evidence.
    soil_layer = row.get("soil_layer") if "soil_layer" in row.index else None
    if pd.notna(soil_layer) and str(soil_layer).strip().lower() == "active layer":
        return (
            "active_layer",
            "active_layer",
            "harmonized_soil_layer",
            f"soil_layer={str(soil_layer).strip()}",
            "strong_inference",
        )

    # Regression/general rule for retrogressive thaw-slump morphology. These
    # labels describe disturbed slump terrain at sample level; generic controls
    # are intentionally left unresolved without publication/site context.
    slump_context = any(term in project_text for term in [
        "retrogressive thaw slump", "retrogressive-thaw-slump", "thaw slump"
    ])
    thaw_feature_terms = [
        "exposed patch", "exposed patches", "drape patch", "drape patches",
        "vegetated raft", "vegetated rafts", "slump floor", "headwall"
    ]
    if slump_context:
        matched = next((term for term in thaw_feature_terms if term in sample_text), None)
        if matched:
            field, value = next(((f, v) for f, v in sample_pieces if matched in v.lower()), (None, None))
            return ("degrading_permafrost", "thaw_affected", "NCBI_BioSample",
                    f"{field}={value}; project_context=retrogressive thaw slump", "strong_inference")

    # A bare occurrence of "permafrost" is not enough to call a sample intact.
    # Likewise, project-level mentions of thermokarst/active layer/permafrost
    # describe study scope, not the state of every sample.
    return "unknown", "exclude", None, None, "unknown"

def score_review_candidate(row):
    """Score scientific usefulness for the convergence-vs-divergence hypothesis.

    The score intentionally weights comparison design more heavily than metadata
    completeness. Maximum nominal score = 30.
    """
    score = 0

    # 1. Within-site comparison coverage (dominant component; max 12).
    coverage = int(row.get("comparison_coverage") or 0)
    if coverage >= 3:
        score += 12
    elif coverage == 2:
        states = set(str(row.get("site_comparison_states") or "").split("|"))
        # Comparisons containing thaw-affected material directly test the primary question.
        score += 8 if "thaw_affected" in states else 5
    elif coverage == 1:
        score += 1

    # 2. Persistent-thaw representation (max 4).
    states = set(str(row.get("site_comparison_states") or "").split("|"))
    if "thaw_affected" in states:
        score += 4

    # 3. Biological replication in this analytical state (max 4).
    nrep = numeric_value(row.get("biological_replicates")) or 0
    if nrep >= 5:
        score += 4
    elif nrep >= 3:
        score += 3
    elif nrep >= 2:
        score += 1

    # 4. Depth comparability among states at the site (max 4).
    depth_comp = row.get("depth_comparability")
    if depth_comp == "good":
        score += 4
    elif depth_comp == "moderate":
        score += 2

    # 5. Technical and provenance quality (max 4).
    if str(row.get("library_source", "")).upper() == "METAGENOMIC": score += 1
    if str(row.get("library_strategy", "")).upper() == "WGS": score += 1
    if row.get("state_confidence") in {"direct", "strong_inference"}: score += 1
    if pd.notna(row.get("latitude")) and pd.notna(row.get("longitude")): score += 1

    # 6. Useful environmental covariates (max 2).
    if pd.notna(row.get("ph")): score += 1
    if pd.notna(row.get("moisture_percent")): score += 1

    return score

def add_hypothesis_design_metrics(review):
    """Annotate within-site state coverage and approximate depth comparability."""
    review = review.copy()
    valid_states = {"intact_permafrost", "active_layer", "thaw_affected"}

    # Site key: prefer submitted site; fall back to rounded coordinates, then project.
    # A row lookup can return a Series when upstream harmonization supplies
    # duplicate column labels. Normalize it before testing missingness.
    def scalar_field(r, field):
        value = r.get(field)
        if isinstance(value, pd.Series):
            nonempty = value.dropna()
            value = nonempty.iloc[0] if not nonempty.empty else None
        return value

    def site_key(r):
        site = scalar_field(r, "site_name")
        if site is not None and pd.notna(site) and str(site).strip():
            return f"site:{str(site).strip().lower()}"
        lat = pd.to_numeric(scalar_field(r, "latitude"), errors="coerce")
        lon = pd.to_numeric(scalar_field(r, "longitude"), errors="coerce")
        if pd.notna(lat) and pd.notna(lon):
            return f"coord:{float(lat):.3f},{float(lon):.3f}"
        return f"project:{scalar_field(r, 'bioproject_accession')}"

    review["_site_key"] = review.apply(site_key, axis=1)
    review["site_comparison_states"] = ""
    review["comparison_coverage"] = 0
    review["depth_comparability"] = "unknown"

    for (_, site), idx in review.groupby(["bioproject_accession", "_site_key"]).groups.items():
        g = review.loc[idx]
        states = sorted(set(g["analysis_state"]).intersection(valid_states))
        state_string = "|".join(states)
        review.loc[idx, "site_comparison_states"] = state_string
        review.loc[idx, "comparison_coverage"] = len(states)

        # Compare median sample mid-depths among represented states. This is a screening
        # metric, not a substitute for study-specific stratigraphic review.
        mids = {}
        for state in states:
            sg = g[g["analysis_state"] == state]
            mid = (pd.to_numeric(sg["depth_cm_min"], errors="coerce") +
                   pd.to_numeric(sg["depth_cm_max"], errors="coerce")) / 2
            if mid.notna().any():
                mids[state] = float(mid.median())
        if len(mids) >= 2:
            spread = max(mids.values()) - min(mids.values())
            if spread <= 20:
                comp = "good"
            elif spread <= 50:
                comp = "moderate"
            else:
                comp = "poor"
            review.loc[idx, "depth_comparability"] = comp

    return review.drop(columns=["_site_key"])

def build_candidate_review(harmonized_df):
    """Create one schema-controlled review row per BioSample x SRA run.

    Ecological eligibility is based on the standardized output of
    eco_harmonization.apply_harmonization(), not raw submitter text alone.
    """
    if harmonized_df.empty:
        return pd.DataFrame(columns=CANDIDATE_REVIEW_COLUMNS)

    records = []
    for _, src in harmonized_df.iterrows():
        rec = {c: None for c in CANDIDATE_REVIEW_COLUMNS}
        rec.update({
            "bioproject_accession": src.get("bioproject_accession"),
            "project_title": src.get("project_title"),
            "biosample_accession": src.get("biosample_accession"),
            "sample_name": first_nonempty(src, ["sample_name", "sample_title"])[0],
            "run_accession": src.get("run_accession"),
            "experiment_accession": src.get("experiment_accession"),
            "library_strategy": src.get("library_strategy"),
            "library_source": src.get("library_source"),
            "platform": src.get("platform"),
            "instrument_model": src.get("instrument_model"),
            "layout": src.get("library_layout"),
            "read_count": numeric_value(src.get("total_spots")),
            "base_count": numeric_value(src.get("total_bases")),
            "site_name": src.get("site_name"),
            "country_region": src.get("location_region"),
            "latitude": src.get("latitude_harmonized"),
            "longitude": src.get("longitude_harmonized"),
            "collection_date": first_nonempty(src, ["collection_date", "sample_collection_date"])[0],
            "ecosystem": src.get("eco_category"),
            "soil_type": src.get("soil_type_reported") or src.get("soil_context"),
            "temperature_c": numeric_value(first_nonempty(src, ["temperature", "temp", "soil_temperature", "sample_temperature"])[0]),
            "ph": numeric_value(first_nonempty(src, ["ph", "soil_ph"])[0]),
            "moisture_percent": numeric_value(first_nonempty(src, ["soil_moisture", "moisture", "water_content"])[0]),
        })

        rec["depth_raw"] = first_nonempty(src, ["depth", "sample_depth", "soil_depth"])[0]
        rec["depth_cm_min"] = src.get("depth_cm_min")
        rec["depth_cm_max"] = src.get("depth_cm_max")
        rec["state_raw"] = first_nonempty(src, ["permafrost_state", "sample_type", "env_local_scale", "env_medium", "isolation_source", "soil_layer", "sample_title", "sample_comment"])[0]
        (rec["permafrost_state"], rec["analysis_state"], rec["state_source"],
         rec["state_evidence"], rec["state_confidence"]) = classify_permafrost_state(src)

        reasons = []
        if str(rec["library_strategy"] or "").upper() != "WGS":
            reasons.append("library_strategy_not_WGS")
        if str(rec["library_source"] or "").upper() != "METAGENOMIC":
            reasons.append("library_source_not_METAGENOMIC")
        if not rec["biosample_accession"]:
            reasons.append("missing_BioSample_accession")

        # Ecological exclusions use harmonized classifications.
        eco_category = str(src.get("eco_category") or "unknown")
        is_soil = src.get("is_soil_sample")
        if is_soil is False or str(is_soil).lower() == "false":
            reasons.append(f"non_soil_ecology:{eco_category}")
        elif eco_category == "unknown":
            reasons.append("ecology_unknown")

        if rec["permafrost_state"] == "unknown":
            reasons.append("permafrost_state_unknown")

        rec["exclusion_reason"] = "; ".join(reasons) if reasons else None

        harmonization_review = any(bool(src.get(flag)) for flag in [
            "eco_review_flag", "location_review_flag", "soil_review_flag"
        ] if pd.notna(src.get(flag)))
        rec["manual_review"] = bool(reasons or rec["analysis_state"] == "exclude" or harmonization_review)
        records.append(rec)

    review = pd.DataFrame(records, columns=CANDIDATE_REVIEW_COLUMNS)

    eligible_state = review["analysis_state"].isin({"intact_permafrost", "active_layer", "thaw_affected"})
    review["biological_replicates"] = float("nan")
    review.loc[eligible_state, "biological_replicates"] = (
        review.loc[eligible_state]
        .groupby(["bioproject_accession", "analysis_state"])["biosample_accession"]
        .transform("nunique")
    )

    review = add_hypothesis_design_metrics(review)

    technical_ok = (
        review["library_strategy"].fillna("").str.upper().eq("WGS")
        & review["library_source"].fillna("").str.upper().eq("METAGENOMIC")
    )
    ecological_ok = ~review["exclusion_reason"].fillna("").str.contains(
        r"non_soil_ecology|ecology_unknown", regex=True
    )
    state_ok = review["analysis_state"].isin({"intact_permafrost", "active_layer", "thaw_affected"})
    repl_ok = review["biological_replicates"].fillna(0).ge(3)

    review["eligibility_tier"] = "C"
    review.loc[technical_ok & ecological_ok & state_ok, "eligibility_tier"] = "B"
    review.loc[technical_ok & ecological_ok & state_ok & repl_ok & review["comparison_coverage"].ge(2), "eligibility_tier"] = "A"
    review.loc[technical_ok & ecological_ok & state_ok & repl_ok & review["comparison_coverage"].ge(3), "eligibility_tier"] = "A+"
    review.loc[~technical_ok | ~ecological_ok, "eligibility_tier"] = "EXCLUDED"
    review["candidate_score"] = review.apply(score_review_candidate, axis=1)

    ordered = CANDIDATE_REVIEW_COLUMNS.copy()
    ordered.insert(ordered.index("candidate_score"), "biological_replicates")
    return review[ordered]


def run_classifier_regression_checks(review):
    """Guard against the two blanket-classification failures found in the first audit."""
    checks = []

    # PRJEB113436 contains multiple terrestrial/aquatic states in project prose,
    # so it must not be blanket-labeled from BioProject text. Harmonized
    # sample-level active-layer calls are allowed.
    x = review[review["bioproject_accession"].eq("PRJEB113436")]
    if not x.empty:
        assert not (x["analysis_state"].isin({"thaw_affected", "intact_permafrost"})).any(), (
            "Regression failure: PRJEB113436 was blanket-classified from project-level text"
        )
        active = x[x["analysis_state"].eq("active_layer")]
        if not active.empty:
            assert active["state_source"].eq("harmonized_soil_layer").all(), (
                "Regression failure: PRJEB113436 active-layer calls must come from harmonized soil-layer evidence"
            )
        checks.append("PRJEB113436: PASS (project prose not propagated; harmonized active layer allowed)")

    # PRJNA1492241 is a retrogressive thaw-slump study. Sample-level disturbed
    # features should be thaw-affected, while generic project 'permafrost' text
    # must not turn controls/other samples into intact permafrost automatically.
    x = review[review["bioproject_accession"].eq("PRJNA1492241")]
    if not x.empty:
        assert (x["analysis_state"] == "thaw_affected").any(), (
            "Regression failure: PRJNA1492241 thaw-slump features were not recognized"
        )
        assert not (x["analysis_state"] == "intact_permafrost").any(), (
            "Regression failure: PRJNA1492241 was labeled intact from project-level permafrost text"
        )
        checks.append("PRJNA1492241: PASS (slump features recognized; no blanket intact label)")

    return checks

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():

    candidate_uids = search_bioprojects(SEARCH_TERM)

    if not candidate_uids:
        print("No candidate BioProjects were returned.")
        return

    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

    # 1. BioProject metadata
    projects_df = retrieve_candidate_metadata(candidate_uids)

    project_columns = [
        "bioproject_accession",
        "bioproject_uid",
        "title",
        "description",
        "organism",
        "taxonomy_id",
        "environment",
        "project_type",
        "submission_date",
        "umbrella_projects",
        "archive_id",
        "retrieval_error",
    ]
    projects_df = projects_df[project_columns]
    projects_df.to_csv(OUTPUT_PROJECTS, index=False)

    print("\nCandidate BioProjects:\n")
    print(
        projects_df[
            ["bioproject_accession", "bioproject_uid", "title", "submission_date"]
        ].to_string(index=False)
    )
    print(f"\nSaved project metadata to: {OUTPUT_PROJECTS}")

    # 2. BioSample metadata
    print("\nRetrieving linked BioSamples...")
    biosamples_df, attributes_df = retrieve_biosamples(projects_df)
    biosamples_df.to_csv(OUTPUT_BIOSAMPLES, index=False)
    attributes_df.to_csv(OUTPUT_ATTRIBUTES, index=False)
    print(f"Saved BioSample table to: {OUTPUT_BIOSAMPLES}")
    print(f"Saved lossless BioSample attribute table to: {OUTPUT_ATTRIBUTES}")

    # 3. SRA experiment/run metadata
    print("\nRetrieving linked SRA records...")
    sra_df = retrieve_sra(projects_df)
    sra_df = filter_sra_records(sra_df)
    sra_df.to_csv(OUTPUT_SRA, index=False)
    print(f"Saved SRA run table to: {OUTPUT_SRA}")

    # 4. Joined CAZyme screening table
    screening_df = build_screening_table(projects_df, biosamples_df, sra_df)
    screening_df.to_csv(OUTPUT_SCREENING, index=False)
    print(f"Saved joined CAZyme screening table to: {OUTPUT_SCREENING}")

    # 5. Ecological/location/soil harmonization
    harmonized_df = apply_harmonization(screening_df)
    harmonized_df.to_csv(OUTPUT_HARMONIZED, index=False)
    print(f"Saved harmonized screening table to: {OUTPUT_HARMONIZED}")

    # 6. Schema-controlled candidate review table
    candidate_review = build_candidate_review(harmonized_df)
    candidate_review.to_csv(OUTPUT_CANDIDATE_REVIEW, index=False)
    print(f"Saved candidate review table to: {OUTPUT_CANDIDATE_REVIEW}")

    print("\nSummary")
    print(f"  BioProjects: {len(projects_df)}")
    print(f"  BioSamples:  {len(biosamples_df)}")
    print(f"  SRA runs:    {len(sra_df)}")

    if not screening_df.empty:
        technical = int(screening_df["cazyme_technical_candidate"].sum())
        print(f"  WGS + METAGENOMIC runs: {technical}")

    if not candidate_review.empty:
        print("  Candidate review tiers:")
        print(candidate_review["eligibility_tier"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
