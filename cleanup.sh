#!/bin/bash

# EcoliTyper Quick Cleanup Script
echo "🧹 EcoliTyper Quick Cleanup"

# Remove build directories
echo "📁 Removing build directories..."
rm -rf build/
rm -rf dist/
rm -rf *.egg-info/
rm -rf ecoliTyper.egg-info/
rm -rf ecolityper.egg-info/

# Remove Python cache files
echo "🐍 Removing Python cache files..."
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find . -type f -name "*.pyc" -delete
find . -type f -name "*.pyo" -delete
find . -type f -name "*.pyd" -delete

# Remove IDE files
echo "💻 Removing IDE files..."
rm -rf .vscode/
rm -rf .idea/
find . -type f -name "*.swp" -delete
find . -type f -name "*.swo" -delete
find . -type f -name "*~" -delete

# Remove test files
echo "🧪 Removing test files..."
rm -f .coverage
rm -rf htmlcov/
rm -rf .pytest_cache/
rm 

echo "✅ Cleanup completed!"
