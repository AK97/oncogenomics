# Source of truth for 96-channel order

BASES = ["A", "C", "G", "T"]
SUBSTITUTIONS = [("C","A"), ("C","G"), ("C","T"), ("T","A"), ("T","C"), ("T","G")]

CHANNELS_96 = [
    f"{left}[{ref}>{alt}]{right}"
    for (ref, alt) in SUBSTITUTIONS
    for left in BASES
    for right in BASES
]