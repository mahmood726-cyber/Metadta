"""Static contracts for the Dose Response Pro artifact layout."""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
APP_HTML = REPO_ROOT / "dose-response-pro.html"
INDEX_HTML = REPO_ROOT / "index.html"
PROTOCOL = REPO_ROOT / "E156-PROTOCOL.md"
RELEASE_NOTES = REPO_ROOT / "release_notes.txt"

EXPECTED_ARTIFACTS = [
    "validate_dose_response_pro.R",
    "r_validation_results.json",
    "r_validation_results.txt",
    "r_simulation_benchmark_results.json",
    "cross_package_benchmark_results.json",
    "strict_beat_r_benchmark_results.json",
    "strict_r_batch_results.json",
    "unit_tests.html",
]


def test_release_notes_and_protocol_align_with_primary_app() -> None:
    html = APP_HTML.read_text(encoding="utf-8")
    protocol = PROTOCOL.read_text(encoding="utf-8")
    release_notes = RELEASE_NOTES.read_text(encoding="utf-8")

    assert "official Dose Response Pro application" in release_notes
    assert "# E156 Protocol" in protocol
    assert "Metadta" in protocol
    assert "<title>Dose Response Pro v18.1" in html
    assert "Dose Response Pro v18.1" in html


def test_external_validation_artifacts_ship_with_repo_and_html_contract() -> None:
    html = APP_HTML.read_text(encoding="utf-8")

    assert "EXTERNAL_ARTIFACT_PATHS" in html
    assert "fetchTextArtifactWithFallback" in html

    expected_root_paths = [
        "validate_dose_response_pro.R",
        "r_validation_results.json",
        "r_validation_results.txt",
        "r_simulation_benchmark_results.json",
    ]
    for path in expected_root_paths:
        assert f"'{path}'" in html or f'"{path}"' in html
        assert (REPO_ROOT / path).is_file(), f"missing artifact: {path}"


def test_core_analysis_shell_and_shipped_artifacts_exist() -> None:
    html = APP_HTML.read_text(encoding="utf-8")
    index_html = INDEX_HTML.read_text(encoding="utf-8")

    expected_markers = [
        "runValidationSuite()",
        "runSimulationBenchmarkGrid()",
        "loadRValidationArtifacts()",
        "loadExternalSimulationBenchmark()",
        "estimateTau2REML(",
        "subgroupForestPlot",
        "forest plot",
        "REML",
    ]
    for marker in expected_markers:
        assert marker in html

    for artifact in EXPECTED_ARTIFACTS:
        assert (REPO_ROOT / artifact).is_file(), f"missing shipped artifact: {artifact}"

    assert "github.com/mahmood726-cyber/Metadta" in index_html
    assert "E156-PROTOCOL.md" in index_html
