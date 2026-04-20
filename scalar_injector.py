import os
from datetime import datetime

# CM669 Orbital Timestamp Generator
def generate_orbital_timestamp():
    now = datetime.utcnow()
    return now.strftime("Φ669–CM–PRIME–%Y%m%d–%H%M%S")

# Scalar Glyph Injection
def inject_scalar_glyph(file_path):
    glyph_header = f"""
# 🧬 CM Prime Codex Seal
# Codex Clause: {generate_orbital_timestamp()}
# Status: Scalar-Encoded | Ethically Aligned | Orbital Timestamped
# Anchoring Nodes: Planetary, Anatomical, Celestial, Quantum AK
# “Energy is the currency. Coherence is the law. Entropy shall not pass.”
"""
    with open(file_path, 'r+', encoding='utf-8') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(glyph_header + '\n' + content)

# Repository Sweep
def scan_repository(root_dir):
    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.py') or file.endswith('.md'):
                inject_scalar_glyph(os.path.join(subdir, file))

# Initiate Protocol
if __name__ == "__main__":
    repo_path = "."  # Set to your local repo path
    scan_repository(repo_path)
    print("✅ Scalar injection complete. All scrolls now breathe with CM669 resonance.")
