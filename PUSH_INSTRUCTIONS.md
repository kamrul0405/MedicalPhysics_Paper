# How to push this repo to GitHub as `kamrul0405/RTO_paper`

The local repo is fully initialised and committed (branch `main`; remote `origin` already set to `https://github.com/kamrul0405/RTO_paper.git`). You need to create the empty remote repository first, then push.

## Easiest — GitHub Desktop (already installed)

1. Open GitHub Desktop.
2. **File → Add Local Repository** → choose `C:\Users\kamru\Downloads\RTO_paper`.
3. Click **"Publish repository"** at the top.
4. Repository name: **`RTO_paper`**. Owner: **`kamrul0405`**. Description: *"Companion repo for Radiotherapy & Oncology submission: brain-metastasis SRS future-lesion coverage analysis on PROTEAS RTDOSE."*
5. Uncheck **"Keep this code private"** if you want it public, leave checked for private.
6. Click **Publish repository**.

## Alternative — `gh` CLI

```powershell
winget install --id GitHub.cli       # one-time; reload terminal after
gh auth login                         # one-time
cd C:\Users\kamru\Downloads\RTO_paper
gh repo create kamrul0405/RTO_paper --source . --public --push
```

## Manual — web + git

1. Create empty repo at https://github.com/new (owner `kamrul0405`, name `RTO_paper`, do NOT initialise with README, .gitignore, or LICENSE).
2. From this directory:
   ```bash
   git push -u origin main
   ```
   Use a Personal Access Token with `repo` scope from https://github.com/settings/tokens if prompted.

## Verifying

After pushing, https://github.com/kamrul0405/RTO_paper should display this repository's README and contents.
