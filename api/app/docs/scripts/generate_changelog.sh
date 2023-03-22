#!/usr/bin/env bash
# PREREQUISITS:
# yarn add -D --modules-folder ./node_modules git-release-notes

DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/../"
GRN="${DIR}/node_modules/.bin/git-release-notes"
CHANGELOG="${DIR}/docs/CHANGELOG.md"

if [ ! -f ${GRN} ]; then
  echo "Installing git-release-notes because it's not there..."
  yarn add -D --modules-folder ${DIR}/node_modules git-release-notes
fi

first=$(git log --reverse --pretty="%h" | head -1)
last=HEAD




echo "Generating CHANGELOG..."
${GRN} -p ${DIR} -b main ${first}..${last} ${DIR}/templates/CHANGELOG.ejs > ${CHANGELOG}


echo "CHANGELOG generated at ${CHANGELOG}."