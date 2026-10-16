#!/usr/bin/env bash
# Usage: ./download_structures.sh <data_dir> train af
# Downloads sequentially with a short pause and retries (the servers answer
# "429 Too Many Requests" to bursts). Existing files are skipped, so it is safe to re-run.
DATA_DIR=$1  # e.g. ../data
SPLIT=$2     # train / eval / test
DB=$3        # pdb / af

OUT="${DATA_DIR}/structures/${SPLIT}/${DB}"
IDS="${DATA_DIR}/splits/${SPLIT}_${DB}.txt"
mkdir -p "$OUT"

download() {  # download <url>; retries with growing backoff
    local url=$1 dest="$OUT/$(basename "$1")"
    [[ -s "$dest" ]] && return 0
    for attempt in 1 2 3 4 5 6; do
        if curl -sf -A "Wget/1.21" -o "$dest.part" "$url"; then
            mv "$dest.part" "$dest"
            sleep 0.2
            return 0
        fi
        sleep $((attempt * 10))
    done
    rm -f "$dest.part"
    echo "FAILED: $url" >&2
    return 1
}

if [[ "$DB" == "pdb" ]]; then
    urls=$(awk -F'_' '{print "https://files.rcsb.org/download/"$1".cif"}' "$IDS" | sort -u)
elif [[ "$DB" == "af" ]]; then
    urls=$(python3 -c "
import re
for line in open('${IDS}'):
    m = re.match(r'AF_AF(.+?)(F\d+)_\d+$', line.strip())
    if m:
        print(f'https://alphafold.ebi.ac.uk/files/AF-{m.group(1)}-{m.group(2)}-model_v6.cif')
" | sort -u)
fi

failed=0
for url in $urls; do
    download "$url" || failed=$((failed + 1))
done
echo "${SPLIT}/${DB}: done, ${failed} failed"
