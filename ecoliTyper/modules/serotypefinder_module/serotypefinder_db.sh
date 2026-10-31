#!/usr/bin/env bash
# serotypefinder_db.sh – download and index the SerotypeFinder database
# Author: Brown Beckley (adapted)
# Usage: ./serotypefinder_db.sh [output_dir]

set -euo pipefail

# Default output directory (relative to script location)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${1:-$SCRIPT_DIR/serotypefinder_db}"

echo "🧬 SerotypeFinder Database Downloader & Indexer"
echo "   Target directory: $OUTDIR"

# Check if Git is available
if ! command -v git &> /dev/null; then
    echo "❌ Error: git is not installed. Please install git first."
    exit 1
fi

# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Error: python3 is not installed. Please install Python 3."
    exit 1
fi

# If the directory exists and contains INSTALL.py, try to update; otherwise clone
if [[ -d "$OUTDIR" ]] && [[ -f "$OUTDIR/INSTALL.py" ]]; then
    echo "📂 Existing SerotypeFinder database found. Updating..."
    cd "$OUTDIR"
    # Pull latest changes (if it's a git repo)
    if [[ -d ".git" ]]; then
        git pull origin master || git pull origin main || echo "⚠️  Git pull failed; continuing with existing files."
    fi
    cd - > /dev/null
else
    if [[ -d "$OUTDIR" ]]; then
        echo "⚠️  Directory '$OUTDIR' exists but INSTALL.py not found. Removing it and cloning fresh..."
        rm -rf "$OUTDIR"
    fi
    echo "📦 Cloning SerotypeFinder database repository..."
    git clone https://bitbucket.org/genomicepidemiology/serotypefinder_db.git "$OUTDIR"
fi

# Check if KMA is installed (required for indexing)
if ! command -v kma &> /dev/null; then
    echo "⚠️  KMA not found. The SerotypeFinder database needs KMA to create indices."
    echo "   Please install KMA (https://bitbucket.org/genomicepidemiology/kma) or ensure it's in your PATH."
    echo "   You can still manually run: cd $OUTDIR && python3 INSTALL.py kma_index"
else
    echo "🔨 Building KMA indices for SerotypeFinder database..."
    cd "$OUTDIR"
    python3 INSTALL.py kma_index
    cd - > /dev/null
    echo "✅ KMA indexing completed."
fi

# Verify important files exist
if [[ -f "$OUTDIR/O_type.fsa" ]] && [[ -f "$OUTDIR/H_type.fsa" ]]; then
    echo "✅ SerotypeFinder database files present."
else
    echo "⚠️  Some expected files (O_type.fsa, H_type.fsa) are missing. The repository structure may have changed."
fi

# Print final status
echo ""
echo "🎉 SerotypeFinder database is ready in: $OUTDIR"
echo ""
echo "📌 Next steps for EcoliTyper:"
echo "   - Set the environment variable: export SEROTYPEFINDER_DB=\"$OUTDIR\""
echo "   - Or update your SerotypeFinder module configuration to point to this directory."
echo ""
echo "🔍 To test the database, run:"
echo "   python3 -c \"import os; print(os.environ.get('SEROTYPEFINDER_DB', 'Not set'))\""
echo ""
echo "✅ Done."
