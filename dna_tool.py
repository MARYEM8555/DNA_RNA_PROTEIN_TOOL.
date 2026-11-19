#!/usr/bin/env python3
"""ADN -> ARNm (transcription) -> Protéine (traduction)"""

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

CODON_TABLE: dict[str, str] = {
    a + b + c: AMINO[i]
    for i, (a, b, c) in enumerate(
        (x, y, z) for x in BASES for y in BASES for z in BASES
    )
}

DNA_COMPLEMENT = str.maketrans("ATCG", "TAGC")


def clean(seq: str) -> str:
    """Remove whitespace/newlines and convert sequence to uppercase."""
    return "".join(seq.upper().split())


def validate_dna(seq: str) -> None:
    """Validate that DNA sequence contains only A, T, C, or G."""
    bad = set(seq) - set("ATCG")
    if bad:
        raise ValueError(f"Invalid DNA bases: {', '.join(sorted(bad))}")


def validate_rna(seq: str) -> None:
    """Validate that RNA sequence contains only A, U, C, or G."""
    bad = set(seq) - set("AUCG")
    if bad:
        raise ValueError(f"Invalid RNA bases: {', '.join(sorted(bad))}")


def reverse_complement(dna: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    cleaned = clean(dna)
    validate_dna(cleaned)
    return cleaned.translate(DNA_COMPLEMENT)[::-1]


def transcribe(dna: str, template_strand: bool = False) -> str:
    """
    Transcribe DNA to mRNA.
    - template_strand=False: coding strand (replace T with U).
    - template_strand=True: template strand (reverse complement, then T -> U).
    """
    cleaned = clean(dna)
    validate_dna(cleaned)
    if template_strand:
        seq = cleaned.translate(DNA_COMPLEMENT)[::-1]
    else:
        seq = cleaned
    return seq.replace("T", "U")


def translate(mrna: str, to_stop: bool = True) -> str:
    """Translate mRNA to protein starting from the first AUG start codon."""
    cleaned = clean(mrna)
    validate_rna(cleaned)

    dna_like = cleaned.replace("U", "T")
    start = dna_like.find("ATG")
    if start == -1:
        raise ValueError("No AUG start codon found in sequence")

    protein: list[str] = []
    for i in range(start, len(dna_like) - 2, 3):
        aa = CODON_TABLE[dna_like[i : i + 3]]
        if aa == "*":
            if to_stop:
                break
            protein.append("*")
        else:
            protein.append(aa)

    return "".join(protein)
