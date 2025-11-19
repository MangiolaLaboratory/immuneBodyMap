#!/bin/bash
# Sync script to update Phoenix_pipeline from HPC_posterior/TAR_SCRIPTS
# Usage: ./sync_phoenix_pipeline.sh

SOURCE_DIR="/home/a1237163/lab/chen/HPC_posterior/TAR_SCRIPTS"
TARGET_DIR="/home/a1237163/lab/chen/immuneHealthyBodyMap/Phoenix_pipeline"

if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory $SOURCE_DIR does not exist"
    exit 1
fi

echo "Syncing files from $SOURCE_DIR to $TARGET_DIR..."

# Remove existing files (but keep .git if it exists)
cd "$TARGET_DIR" || exit 1
find . -type f ! -path './.git/*' -delete
find . -type d -empty ! -path './.git/*' -delete

# Copy all files from source
cp -r "$SOURCE_DIR"/* "$TARGET_DIR"/ 2>/dev/null
cp -r "$SOURCE_DIR"/. "$TARGET_DIR"/ 2>/dev/null 2>&1 | grep -v "cannot stat"

echo "Sync complete!"
echo "Files in Phoenix_pipeline have been updated from HPC_posterior/TAR_SCRIPTS"
echo ""
echo "You can now:"
echo "  1. Review changes: git diff Phoenix_pipeline/"
echo "  2. Stage changes: git add Phoenix_pipeline/"
echo "  3. Commit: git commit -m 'Update Phoenix_pipeline from source'"
echo "  4. Push: git push"

