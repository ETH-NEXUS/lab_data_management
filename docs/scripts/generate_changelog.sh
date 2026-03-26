#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)/../"
CHANGELOG="${DIR}/docs/CHANGELOG.md"

if ! command -v git-release-notes >/dev/null 2>&1; then
  echo "git-release-notes is not installed in the docs container."
  exit 1
fi

first=$(git -C "${DIR}" log --reverse --pretty="%h" | head -1)
last=HEAD

echo "Generating CHANGELOG..."
git-release-notes -p "${DIR}" -b main "${first}..${last}" "${DIR}/templates/CHANGELOG.ejs" > "${CHANGELOG}"

echo "CHANGELOG generated at ${CHANGELOG}."
