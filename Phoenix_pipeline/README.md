# Phoenix Pipeline

This directory contains the pipeline scripts for Phoenix HPC analysis.

## Syncing with Source Directory

This directory is synced from `/home/a1237163/lab/chen/HPC_posterior/TAR_SCRIPTS`.

To update this directory after making changes to the source:

```bash
# From the immuneHealthyBodyMap directory
./sync_phoenix_pipeline.sh
```

Then commit and push the changes:

```bash
git add Phoenix_pipeline/
git commit -m "Update Phoenix_pipeline from source"
git push
```

## Note

The files in this directory are actual copies (not symlinks) so they can be properly tracked by Git and pushed to GitHub. Use the sync script to keep them updated with the source directory.

