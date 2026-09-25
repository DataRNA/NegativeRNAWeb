from pathlib import Path

from docx import Document
from docx.shared import RGBColor


SOURCE = Path(
    r"C:\Users\DataRNA1\Desktop\makale\nerna_article_BMC revised_final_with_reviewer_comments.docx"
)
OUTPUT = Path(
    r"C:\Users\DataRNA1\Desktop\makale\nerna_article_BMC revised_final_without_motif_metric.docx"
)

REPLACEMENTS = {
    (
        "Therefore, negative data generation for ncRNA machine learning is not simply a "
        "technical convenience; it is a core methodological challenge that affects the model "
        "calibration and biological interpretability. Simple random shuffling may distort "
        "higher-order dependencies, whereas dinucleotide-preserving shuffling maintains local "
        "composition but may still produce unrealistic structural profiles or disrupt motifs "
        "in ways that do not reflect plausible biological variation. Conversely, negatives "
        "constructed from unrelated genomic regions may differ systematically in length, GC "
        "content, or repeat composition, producing models that learn dataset-specific shortcuts "
        "rather than the discriminative features of the ncRNA class. A standardised, "
        "reproducible, and configurable negative data generation strategy is required to enable "
        "fair benchmarking and support robust model development across multiple ncRNA biotypes."
    ): (
        "Therefore, negative data generation for ncRNA machine learning is not simply a "
        "technical convenience; it is a core methodological challenge that affects the model "
        "calibration and biological interpretability. Simple random shuffling may distort "
        "higher-order dependencies, whereas dinucleotide-preserving shuffling maintains local "
        "composition but may still produce unrealistic structural profiles or disrupt "
        "higher-order sequence patterns in ways that do not reflect plausible biological "
        "variation. Conversely, negatives constructed from unrelated genomic regions may differ "
        "systematically in length, GC content, or repeat composition, producing models that "
        "learn dataset-specific shortcuts rather than the discriminative features of the ncRNA "
        "class. A standardised, reproducible, and configurable negative data generation strategy "
        "is required to enable fair benchmarking and support robust model development across "
        "multiple ncRNA biotypes."
    ),
    (
        "Dinucleotide-preserving shuffling was applied to maintain the dinucleotide composition "
        "of each sequence while disrupting the higher-order motifs.\u00a0 Because the dinucleotide "
        "composition is preserved by design, gross compositional properties such as GC content "
        "are expected to remain similar between the original and shuffled sequences, whereas "
        "local motifs and longer patterns may be disrupted."
    ): (
        "Dinucleotide-preserving shuffling was applied to maintain the dinucleotide composition "
        "of each sequence while altering higher-order sequence patterns. Because dinucleotide "
        "composition is preserved by design, gross compositional properties such as GC content "
        "are expected to remain similar between the original and shuffled sequences, whereas "
        "longer sequence patterns may change."
    ),
}


def replace_paragraph_text(paragraph, replacement: str) -> None:
    """Replace one single-run paragraph and mark the revised paragraph in red."""
    if len(paragraph.runs) != 1:
        raise RuntimeError(
            "Expected a single-run paragraph; refusing to rewrite complex Word markup."
        )
    run = paragraph.runs[0]
    run.text = replacement
    run.font.color.rgb = RGBColor(255, 0, 0)


def main() -> None:
    document = Document(SOURCE)
    pending = dict(REPLACEMENTS)

    for paragraph in document.paragraphs:
        replacement = pending.pop(paragraph.text, None)
        if replacement is not None:
            replace_paragraph_text(paragraph, replacement)

    if pending:
        missing = "\n".join(pending)
        raise RuntimeError(f"Target paragraphs were not found:\n{missing}")

    document.save(OUTPUT)
    print(OUTPUT)


if __name__ == "__main__":
    main()
