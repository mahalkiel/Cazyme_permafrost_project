#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
First-pass ecological, location, and soil harmonization for the CAZyme
screening table.

The script keeps all raw screening columns intact and appends derived columns
with confidence/evidence fields. The goal is an auditable review table, not a
final ontology-grade annotation.

Created on Thu Sep 17 09:46:36 2026

@author: michaelsnyder, codex assisted
"""

from pathlib import Path
import re

import pandas as pd


OUTPUT_FOLDER = Path(__file__).resolve().parent.parent / "output" / "explr_metadata"

OUTPUT_SCREENING = OUTPUT_FOLDER / "candidate_cazyme_screening.csv"
OUTPUT_HARMONIZED = OUTPUT_FOLDER / "candidate_cazyme_screening_harmonized.csv"


TEXT_FIELDS = [
    "gold_ecosystem_classification",
    "env_broad_scale",
    "env_local_scale",
    "env_medium",
    "organism",
    "project_organism",
    "package",
    "isolation_source",
    "sample_title",
    "experiment_title",
    "project_title",
    "project_description",
    "sample_comment",
    "geo_loc_name",
]


def clean_text(value):
    """Return a stripped string, or an empty string for missing values."""
    if pd.isna(value):
        return ""
    return str(value).strip()


def row_text(row, fields=TEXT_FIELDS):
    """Combine selected row fields into lowercase searchable text."""
    return " | ".join(clean_text(row.get(field)) for field in fields).lower()


def first_present(row, fields):
    """Return the first non-empty value and its source field name."""
    for field in fields:
        value = clean_text(row.get(field))
        if value:
            return field, value
    return "", ""


def evidence(field, value):
    """Format a compact evidence string for derived labels."""
    if not field or not value:
        return ""
    return f"{field}={value}"


def contains_any(text, terms):
    return any(term in text for term in terms)


def harmonize_ecology(row):
    """
    Classify broad ecological context using structured environmental fields
    first, then less controlled submitter text.
    """
    gold = clean_text(row.get("gold_ecosystem_classification"))
    env_medium = clean_text(row.get("env_medium"))
    env_broad = clean_text(row.get("env_broad_scale"))
    env_local = clean_text(row.get("env_local_scale"))
    isolation = clean_text(row.get("isolation_source"))
    organism = clean_text(row.get("organism")) or clean_text(row.get("project_organism"))
    text = row_text(row)

    result = {
        "eco_domain": "unknown",
        "eco_category": "unknown",
        "eco_subcategory": "",
        "eco_confidence": "low",
        "eco_evidence": "",
        "eco_review_flag": True,
    }

    gold_lower = gold.lower()
    if "aquatic | marine" in gold_lower or "marine | oceanic" in gold_lower:
        result.update({
            "eco_domain": "aquatic",
            "eco_category": "marine",
            "eco_subcategory": "oceanic",
            "eco_confidence": "high",
            "eco_evidence": evidence("gold_ecosystem_classification", gold),
            "eco_review_flag": False,
        })
        return result

    medium_lower = env_medium.lower()
    broad_lower = env_broad.lower()
    local_lower = env_local.lower()
    isolation_lower = isolation.lower()
    organism_lower = organism.lower()

    if "soil" in medium_lower:
        result.update({
            "eco_domain": "terrestrial",
            "eco_category": "soil",
            "eco_subcategory": infer_soil_context(row),
            "eco_confidence": "high",
            "eco_evidence": evidence("env_medium", env_medium),
            "eco_review_flag": False,
        })
        return result

    if "snow" in medium_lower or "snow microcosm" in broad_lower:
        result.update({
            "eco_domain": "cryosphere",
            "eco_category": "snow",
            "eco_subcategory": "snow microcosm" if "microcosm" in text else "",
            "eco_confidence": "high",
            "eco_evidence": evidence("env_medium", env_medium) or evidence("env_broad_scale", env_broad),
            "eco_review_flag": False,
        })
        return result

    if "ice" in medium_lower or "|ice" in local_lower:
        result.update({
            "eco_domain": "cryosphere",
            "eco_category": "ice",
            "eco_subcategory": "water-body-derived ice" if "water-body-derived" in medium_lower else "",
            "eco_confidence": "high",
            "eco_evidence": evidence("env_medium", env_medium) or evidence("env_local_scale", env_local),
            "eco_review_flag": False,
        })
        return result

    if "sediment" in medium_lower or "sediment" in local_lower or "sediment" in isolation_lower:
        result.update({
            "eco_domain": "aquatic",
            "eco_category": "sediment",
            "eco_subcategory": "surface sediment" if "surface" in text else "",
            "eco_confidence": "high" if "sediment" in medium_lower or "sediment" in local_lower else "medium",
            "eco_evidence": (
                evidence("env_medium", env_medium)
                or evidence("env_local_scale", env_local)
                or evidence("isolation_source", isolation)
            ),
            "eco_review_flag": False,
        })
        return result

    if "brackish water" in medium_lower:
        result.update({
            "eco_domain": "aquatic",
            "eco_category": "brackish water",
            "eco_subcategory": "",
            "eco_confidence": "high",
            "eco_evidence": evidence("env_medium", env_medium),
            "eco_review_flag": False,
        })
        return result

    if "pingo biome" in broad_lower:
        result.update({
            "eco_domain": "terrestrial",
            "eco_category": "soil",
            "eco_subcategory": "pingo",
            "eco_confidence": "medium",
            "eco_evidence": evidence("env_broad_scale", env_broad),
            "eco_review_flag": False,
        })
        return result

    if "soil metagenome" in organism_lower:
        result.update({
            "eco_domain": "terrestrial",
            "eco_category": "soil",
            "eco_subcategory": infer_soil_context(row),
            "eco_confidence": "medium",
            "eco_evidence": evidence("organism", organism),
            "eco_review_flag": False,
        })
        return result

    if contains_any(isolation_lower, ["top soil", "subsurface soil", "tundra soil", "surface soil", "arctic tundra soil"]):
        result.update({
            "eco_domain": "terrestrial",
            "eco_category": "soil",
            "eco_subcategory": infer_soil_context(row),
            "eco_confidence": "medium",
            "eco_evidence": evidence("isolation_source", isolation),
            "eco_review_flag": False,
        })
        return result

    if "seawater" in isolation_lower or "seawater" in organism_lower:
        result.update({
            "eco_domain": "aquatic",
            "eco_category": "marine",
            "eco_subcategory": "seawater",
            "eco_confidence": "medium",
            "eco_evidence": evidence("isolation_source", isolation) or evidence("organism", organism),
            "eco_review_flag": False,
        })
        return result

    if "snow" in isolation_lower or "snow metagenome" in organism_lower:
        result.update({
            "eco_domain": "cryosphere",
            "eco_category": "snow",
            "eco_subcategory": "",
            "eco_confidence": "medium",
            "eco_evidence": evidence("isolation_source", isolation) or evidence("organism", organism),
            "eco_review_flag": False,
        })
        return result

    if "ice sample" in isolation_lower or "arctic ice" in isolation_lower:
        result.update({
            "eco_domain": "cryosphere",
            "eco_category": "ice",
            "eco_subcategory": "ice sample",
            "eco_confidence": "medium",
            "eco_evidence": evidence("isolation_source", isolation),
            "eco_review_flag": False,
        })
        return result

    if re.search(r"\bice\d*\b|\bsea-ice\b|\bsea ice\b", text):
        result.update({
            "eco_domain": "cryosphere",
            "eco_category": "ice",
            "eco_subcategory": "sea ice" if "sea ice" in text or "sea-ice" in text else "ice sample",
            "eco_confidence": "medium",
            "eco_evidence": "keyword match in title/description",
            "eco_review_flag": False,
        })
        return result

    if "water sample" in isolation_lower or "water metagenome" in organism_lower:
        result.update({
            "eco_domain": "aquatic",
            "eco_category": "water",
            "eco_subcategory": "arctic water" if "arctic" in text else "",
            "eco_confidence": "medium",
            "eco_evidence": evidence("isolation_source", isolation) or evidence("organism", organism),
            "eco_review_flag": False,
        })
        return result

    source_text = " | ".join(
        clean_text(row.get(field))
        for field in ["isolation_source", "sample_title", "experiment_title", "organism"]
    ).lower()

    if re.search(r"\bwater\b", source_text):
        result.update({
            "eco_domain": "aquatic",
            "eco_category": "water",
            "eco_subcategory": "arctic water" if "arctic" in text else "",
            "eco_confidence": "medium",
            "eco_evidence": "keyword match in source/sample fields",
            "eco_review_flag": False,
        })
        return result

    if "cold spring" in organism_lower:
        result.update({
            "eco_domain": "aquatic",
            "eco_category": "spring",
            "eco_subcategory": "cold spring",
            "eco_confidence": "medium",
            "eco_evidence": evidence("organism", organism),
            "eco_review_flag": False,
        })
        return result

    if contains_any(text, ["soil", "tundra", "gelisol", "histosol", "permafrost"]):
        result.update({
            "eco_domain": "terrestrial",
            "eco_category": "soil",
            "eco_subcategory": infer_soil_context(row),
            "eco_confidence": "low",
            "eco_evidence": "keyword match in title/description",
            "eco_review_flag": True,
        })
        return result

    return result


def normalize_geo_name(value):
    """Clean common location encoding and separator issues without hiding raw data."""
    value = clean_text(value)
    replacements = {
        "Br_ggerbreen": "Broggerbreen",
        "BrÃ¸ggerbreen": "Broggerbreen",
        "Ny-Ã\x85lesund": "Ny-Alesund",
        "Ny Alesund": "Ny-Alesund",
    }
    for old, new in replacements.items():
        value = value.replace(old, new)
    return value


def parse_lat_lon_text(value):
    """
    Parse common BioSample coordinate strings such as
    '78.91517 N 11.82717 E'.
    """
    value = clean_text(value)
    if not value:
        return None, None

    pattern = re.compile(
        r"(?P<lat>-?\d+(?:\.\d+)?)\s*(?P<lat_dir>[NS])?[,\s]+"
        r"(?P<lon>-?\d+(?:\.\d+)?)\s*(?P<lon_dir>[EW])?",
        re.IGNORECASE,
    )
    match = pattern.search(value)
    if not match:
        return None, None

    lat = float(match.group("lat"))
    lon = float(match.group("lon"))
    lat_dir = (match.group("lat_dir") or "").upper()
    lon_dir = (match.group("lon_dir") or "").upper()

    if lat_dir == "S":
        lat = -abs(lat)
    elif lat_dir == "N":
        lat = abs(lat)

    if lon_dir == "W":
        lon = -abs(lon)
    elif lon_dir == "E":
        lon = abs(lon)

    return lat, lon


def first_coordinate(row, canonical_field, alternate_fields):
    """Find the first coordinate value from preferred and alternate columns."""
    for field in [canonical_field] + alternate_fields:
        value = clean_text(row.get(field))
        if not value:
            continue
        try:
            return float(value), field
        except ValueError:
            continue
    return None, ""


def infer_polar_region(geo_name, latitude, text):
    geo_lower = geo_name.lower()
    if "antarctica" in geo_lower or "king george island" in geo_lower:
        return "Antarctic"
    if contains_any(geo_lower, ["svalbard", "arctic ocean", "victoria island", "cambridge bay"]):
        return "Arctic"
    if "tibetan plateau" in geo_lower or "beiluhe" in geo_lower:
        return "alpine permafrost"
    if latitude is not None:
        if latitude <= -60:
            return "Antarctic"
        if latitude >= 66.5:
            return "Arctic"
        if 50 <= latitude < 66.5:
            return "subarctic/boreal"
    if contains_any(text, ["arctic", "svalbard"]):
        return "Arctic"
    if "antarctic" in text:
        return "Antarctic"
    if contains_any(text, ["permafrost", "tibetan plateau"]):
        return "permafrost-associated"
    return "unknown"


def harmonize_location(row):
    geo_name = normalize_geo_name(row.get("geo_loc_name"))
    text = row_text(row)

    lat, lon = parse_lat_lon_text(row.get("lat_lon"))
    coord_evidence = "lat_lon" if lat is not None and lon is not None else ""

    if lat is None:
        lat, lat_field = first_coordinate(row, "latitude", ["geographic_location_latitude"])
    else:
        lat_field = coord_evidence

    if lon is None:
        lon, lon_field = first_coordinate(row, "longitude", ["geographic_location_longitude"])
    else:
        lon_field = coord_evidence

    location_region = "unknown"
    site_name = ""
    if geo_name:
        if ":" in geo_name:
            location_region, site_name = [part.strip() for part in geo_name.split(":", 1)]
        else:
            location_region = geo_name

    hemisphere = "unknown"
    if lat is not None:
        hemisphere = "northern" if lat >= 0 else "southern"

    polar_region = infer_polar_region(geo_name, lat, text)

    confidence = "low"
    if lat is not None and lon is not None:
        confidence = "high"
    elif geo_name:
        confidence = "medium"

    evidence_parts = []
    if geo_name:
        evidence_parts.append(evidence("geo_loc_name", geo_name))
    if coord_evidence:
        evidence_parts.append(evidence("lat_lon", clean_text(row.get("lat_lon"))))
    else:
        if lat_field:
            evidence_parts.append(evidence(lat_field, lat))
        if lon_field:
            evidence_parts.append(evidence(lon_field, lon))

    return {
        "location_region": location_region,
        "site_name": site_name,
        "latitude_harmonized": lat,
        "longitude_harmonized": lon,
        "hemisphere": hemisphere,
        "polar_region": polar_region,
        "location_confidence": confidence,
        "location_evidence": "; ".join(evidence_parts),
        "location_review_flag": confidence == "low" or polar_region == "unknown",
    }


def parse_depth_cm(value):
    """
    Parse common depth strings into centimeter min/max values.

    Returns (min_cm, max_cm). For point values, min and max are the same.
    """
    value = clean_text(value).lower().replace("_", "-")
    if not value or value in {"surface", "missing: control sample"}:
        if value == "surface":
            return 0.0, 0.0
        return None, None

    meter_match = re.fullmatch(r"(\d+(?:\.\d+)?)\s*m", value)
    if meter_match:
        cm = float(meter_match.group(1)) * 100
        return cm, cm

    range_match = re.search(
        r"(?P<min>\d+(?:\.\d+)?)\s*-\s*(?P<max>\d+(?:\.\d+)?)\s*(?P<unit>cm|m)?",
        value,
    )
    if range_match:
        min_depth = float(range_match.group("min"))
        max_depth = float(range_match.group("max"))
        unit = range_match.group("unit") or "cm"
        if unit == "m":
            min_depth *= 100
            max_depth *= 100
        return min_depth, max_depth

    point_match = re.search(r"(?P<value>\d+(?:\.\d+)?)\s*(?P<unit>cm|m)?", value)
    if point_match:
        depth = float(point_match.group("value"))
        if (point_match.group("unit") or "cm") == "m":
            depth *= 100
        return depth, depth

    return None, None


def infer_soil_layer(row):
    text = row_text(row)
    isolation = clean_text(row.get("isolation_source")).lower()
    depth_min, depth_max = parse_depth_cm(row.get("depth"))
    if depth_min is None and depth_max is None:
        depth_min, depth_max = parse_depth_cm(row.get("isolation_source"))

    if "top soil" in isolation or "topsoil" in text:
        return "topsoil"
    if "subsurface soil" in isolation:
        return "subsurface"
    if "surface soil" in isolation or "surface tundra soil" in text:
        return "surface"
    if "active layer" in text:
        return "active layer"
    if depth_min is not None and depth_max is not None:
        if depth_max <= 10:
            return "surface"
        return "subsurface"
    if "permafrost" in text:
        return "permafrost-associated"
    return "unknown"


def infer_soil_type_reported(row):
    text = row_text(row)
    if "gelisol" in text:
        return "Gelisol"
    if "histosol" in text:
        return "Histosol"
    if "cryosol" in text:
        return "Cryosol"
    if "tundra soil" in text:
        return "tundra soil"
    if "arctic surface soil" in text:
        return "Arctic surface soil"
    return ""


def infer_soil_context(row):
    text = row_text(row)
    if "glacier foreland" in text or "austre" in text or "broggerbreen" in text:
        return "glacier foreland"
    if "pingo" in text:
        return "pingo"
    if "tundra with shrubs" in text or "shrub" in text:
        return "shrub tundra"
    if "moss" in text or "sedge" in text or "sphag" in text:
        return "moss/sedge tundra"
    if "tundra" in text:
        return "tundra"
    if "vegetated rafts" in text:
        return "vegetated raft"
    if "exposed patches" in text:
        return "exposed patch"
    if "drape patches" in text:
        return "drape patch"
    if "control patches" in text:
        return "control patch"
    if "water-affected" in text:
        return "water-affected area"
    return ""


def harmonize_soil(row):
    text = row_text(row)
    eco = harmonize_ecology(row)
    depth_min, depth_max = parse_depth_cm(row.get("depth"))
    if depth_min is None and depth_max is None:
        depth_min, depth_max = parse_depth_cm(row.get("isolation_source"))
    soil_type = infer_soil_type_reported(row)
    soil_context = infer_soil_context(row)

    has_soil_terms = contains_any(
        text,
        ["soil", "gelisol", "histosol", "cryosol", "pingo"],
    )
    is_soil = eco["eco_category"] == "soil" or (
        has_soil_terms and eco["eco_category"] not in {"water", "ice", "snow", "marine", "sediment", "brackish water", "spring"}
    )

    confidence = "low"
    evidence_parts = []
    field, value = first_present(row, ["env_medium", "isolation_source", "depth", "project_description"])
    if field:
        evidence_parts.append(evidence(field, value))

    if clean_text(row.get("env_medium")).lower().find("soil") >= 0:
        confidence = "high"
    elif is_soil and clean_text(row.get("isolation_source")):
        confidence = "medium"
    elif is_soil:
        confidence = "low"

    soil_layer = infer_soil_layer(row) if is_soil else ""

    return {
        "is_soil_sample": bool(is_soil),
        "soil_layer": soil_layer,
        "soil_type_reported": soil_type if is_soil else "",
        "soil_context": soil_context if is_soil else "",
        "depth_cm_min": depth_min if is_soil else None,
        "depth_cm_max": depth_max if is_soil else None,
        "soil_confidence": confidence if is_soil else "",
        "soil_evidence": "; ".join(evidence_parts) if is_soil else "",
        "soil_review_flag": bool(is_soil and (confidence == "low" or soil_layer == "unknown")),
    }


def apply_harmonization(screen_df):
    """Append harmonization columns to the screening dataframe."""
    ecology_df = pd.DataFrame([harmonize_ecology(row) for _, row in screen_df.iterrows()])
    location_df = pd.DataFrame([harmonize_location(row) for _, row in screen_df.iterrows()])
    soil_df = pd.DataFrame([harmonize_soil(row) for _, row in screen_df.iterrows()])

    return pd.concat(
        [
            screen_df.reset_index(drop=True),
            ecology_df.reset_index(drop=True),
            location_df.reset_index(drop=True),
            soil_df.reset_index(drop=True),
        ],
        axis=1,
    )


def main():
    print(f"Loading screening table: {OUTPUT_SCREENING}")
    screen_df = pd.read_csv(OUTPUT_SCREENING)

    harmonized_df = apply_harmonization(screen_df)
    harmonized_df.to_csv(OUTPUT_HARMONIZED, index=False)

    print(f"Wrote harmonized table: {OUTPUT_HARMONIZED}")
    print(f"Rows: {len(harmonized_df)}")

    for column in ["eco_category", "polar_region", "soil_layer"]:
        print(f"\n{column}")
        print(harmonized_df[column].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
