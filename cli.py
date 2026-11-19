#!/usr/bin/env python3
"""Command-line interface for DNA / mRNA / Protein tool."""

from dna_tool import reverse_complement, transcribe, translate


def main() -> None:
    print("=== DNA Tool: Transcription & Translation ===")
    while True:
        print("\n1) Transcribe (DNA -> mRNA)")
        print("2) Translate (mRNA -> Protein)")
        print("3) Both (DNA -> mRNA -> Protein)")
        print("4) Reverse complement")
        print("0) Exit")
        choice = input("Choose an option: ").strip()

        try:
            if choice == "0":
                break
            elif choice == "1":
                seq = input("DNA: ")
                tpl = input("Is this a template strand? (y/n): ").strip().lower() == "y"
                print("mRNA:", transcribe(seq, tpl))
            elif choice == "2":
                print("Protein:", translate(input("mRNA: ")))
            elif choice == "3":
                seq = input("DNA: ")
                tpl = input("Is this a template strand? (y/n): ").strip().lower() == "y"
                mrna = transcribe(seq, tpl)
                print("mRNA:", mrna)
                print("Protein:", translate(mrna))
            elif choice == "4":
                print(reverse_complement(input("DNA: ")))
            else:
                print("Invalid choice")
        except ValueError as e:
            print("Error:", e)


if __name__ == "__main__":
    main()
