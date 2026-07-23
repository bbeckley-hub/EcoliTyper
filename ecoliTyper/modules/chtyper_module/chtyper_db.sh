#!/usr/bin/env bash
# chTyper_db.sh – download and prepare the CHTyper database for EcoliTyper
# Author: Brown Beckley (adapted)
# Usage: ./chtyper_db.sh [output_dir]

set -euo pipefail

# Default output directory
OUTDIR="${1:-./chtyper_db}"

echo "🔽 CHTyper Database Downloader"
echo "   Target directory: $OUTDIR"

# Check if Git is available
if ! command -v git &> /dev/null; then
    echo "❌ Error: git is not installed. Please install git first."
    exit 1
fi

# If the directory exists, try to update; otherwise clone
if [[ -d "$OUTDIR/.git" ]]; then
    echo "📂 Existing repository found. Updating..."
    cd "$OUTDIR"
    git pull origin master || git pull origin main
    cd - > /dev/null
else
    if [[ -d "$OUTDIR" ]]; then
        echo "⚠️  Directory '$OUTDIR' exists but is not a git repository."
        echo "   Removing it and cloning fresh..."
        rm -rf "$OUTDIR"
    fi
    echo "📦 Cloning CHTyper database repository..."
    git clone https://bitbucket.org/genomicepidemiology/chtyper_db.git "$OUTDIR"
fi

# Now check if the database files are present
if [[ ! -f "$OUTDIR/chtyper.fasta" ]] && [[ ! -f "$OUTDIR/chtyper_db.fasta" ]]; then
    echo "⚠️  Warning: No FASTA file found in the database directory."
    echo "   The repository structure may have changed. Listing contents:"
    ls -la "$OUTDIR"
else
    echo "✅ CHTyper database downloaded successfully."
fi

# Optional: Build BLAST databases if BLAST+ is available
if command -v makeblastdb &> /dev/null; then
    echo "🧬 Building BLAST databases (if needed)..."
    for fasta in "$OUTDIR"/*.fasta; do
        if [[ -f "$fasta" ]]; then
            base=$(basename "$fasta" .fasta)
            if [[ ! -f "$OUTDIR/${base}.phr" ]]; then
                makeblastdb -in "$fasta" -dbtype nucl -out "$OUTDIR/$base"
                echo "   Built BLAST db for $base"
            fi
        fi
    done
else
    echo "ℹ️  BLAST+ not found; skipping BLAST database build."
    echo "   If your CHTyper module requires BLAST, install BLAST+ and rerun this script."
fi

echo "🎉 Done. The CHTyper database is ready in: $OUTDIR"
echo "   You may need to update your CHTyper configuration to point to this directory."
