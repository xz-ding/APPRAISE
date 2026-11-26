"""Restore chain boundaries in single-chain ESMFold v1 PDB exports."""

from string import ascii_lowercase, ascii_uppercase


_AMINO_ACID = dict(zip(
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR UNK".split(),
    "ACDEFGHIKLMNPQRSTVWYX",
))
_CHAIN_IDS = ascii_uppercase + ascii_lowercase


def linker_length(joined_sequence, chain_sequences):
    """Validate the joined sequence and return its uniform artificial linker length."""
    if len(chain_sequences) < 2 or not all(chain_sequences):
        raise ValueError("Supply the original separated chain sequences in modeling order.")
    extra = len(joined_sequence) - sum(map(len, chain_sequences))
    length, remainder = divmod(extra, len(chain_sequences) - 1)
    if length < 0 or remainder or ("G" * length).join(chain_sequences) != joined_sequence:
        raise ValueError("The joined sequence does not match the original chains separated by glycine linkers.")
    return length


def restore_esmfold_v1_chains(pdb_text, chain_sequences):
    """Use source sequences to split one complete joined ESMFold v1 PDB model.

    The sequences must come from the separated input, in the order used for
    prediction. The input must contain one chain with consecutive residue IDs.
    Atom coordinates and confidence values are preserved.
    """
    if not isinstance(chain_sequences, (list, tuple)) or not 2 <= len(chain_sequences) <= len(_CHAIN_IDS):
        raise ValueError("Supply the original separated chain sequences in modeling order.")
    chains = tuple(sequence.strip().upper() for sequence in chain_sequences)
    lines = pdb_text.splitlines()
    atoms = [line for line in lines if line.startswith("ATOM  ")]
    if (not atoms or len({line[21:22] for line in atoms}) != 1
            or sum(line.startswith("MODEL ") for line in lines) > 1
            or any(len(line) < 66 or line[26] != " " for line in atoms)
            or any(line.startswith("HETATM") for line in lines)):
        raise ValueError("Expected a single joined ESMFold v1 PDB with standard protein residue numbering.")

    residues = {}
    for line in atoms:
        resi = int(line[22:26])
        name = line[17:20].strip()
        if residues.setdefault(resi, name) != name:
            raise ValueError("The PDB contains conflicting residue identities.")
    first = min(residues)
    if list(residues) != list(range(first, first + len(residues))):
        raise ValueError("The PDB residue numbers must be consecutive and in order.")
    try:
        sequence = "".join(_AMINO_ACID[name] for name in residues.values())
    except KeyError as exc:
        raise ValueError(f"Unrecognized ESMFold v1 residue name {exc.args[0]!r}.") from exc
    gap = linker_length(sequence, chains)

    retained = {}
    for chain_id, chain in zip(_CHAIN_IDS, chains):
        for i in range(len(chain)):
            retained[first + i] = (chain_id, i + 1)
        first += len(chain) + gap

    output = []
    previous_chain = None
    for line in atoms:
        layout = retained.get(int(line[22:26]))
        if layout is None:
            continue
        chain_id, resi = layout
        if previous_chain and previous_chain != chain_id:
            output.append("TER")
        output.append(line[:21] + f"{chain_id}{resi:4d}" + line[26:])
        previous_chain = chain_id
    return "\n".join([*output, "TER", "END", ""])
