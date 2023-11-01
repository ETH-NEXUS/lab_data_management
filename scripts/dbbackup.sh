#!/bin/bash


DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LOG_FILE="${DIR}/backup.log"

function log {
    echo "$(date +'%Y-%m-%d %H:%M:%S') $@" | tee -a "$LOG_FILE"
}

log "Starting DB Backup..."
cd "${DIR}/.." && ./manage.py dbbackup
log "DONE DB Backup."
