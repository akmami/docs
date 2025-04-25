# Clean up the Github Repository

When your git repository becomes large in the sense that large files exist in your history that you don't use, installation is a pain in the a...

To overcome this problem, you can purge and re-write the git history.

Here is a script to see all the files created in your repository.

```bash
git rev-list --objects --all | \
    git cat-file --batch-check='%(objecttype) %(objectname) %(rest)' | \
    grep 'blob' | \
    while read type hash filename; do \
        size=$(git cat-file -s "$hash"); \
        size_mb=$(echo "scale=2; $size / 1024 / 1024" | bc); \
        echo "$filename $size_mb MB"; \
    done | sort -k2 -n
```

Now, if you decide to delete some files from history, you can run the following script using the `git-filter-repo` Python package.

You might need to give the binary location of Python language before calling `python3 git-filter-repo ...`

```bash
# remove file
python3 git-filter-repo --path <path-to-large-file> --invert-paths [--force]

# force push the changes
git remote add origin <your-repo-url>
git push --set-upstream origin main [--force]
git push origin --force --all
git push origin --force --tags

# notify people
git reflog expire --expire=now --all
git gc --prune=now --aggressive

# check repo size
git count-objects -vH
```

BAAAAM!!! Nothing better than (almost) a fresh start...
