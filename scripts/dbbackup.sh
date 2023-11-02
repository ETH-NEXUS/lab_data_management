#!/bin/bash


DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LOG_FILE="${DIR}/backup.log"

function log {
    echo "$(date +'%Y-%m-%d %H:%M:%S') $@" | tee -a "$LOG_FILE"
}

log "Starting DB Backup..."
cd "${DIR}/.." && docker-compose exec api /bin/bash -c './manage.py dbbackup' 2>&1
log "DONE DB Backup."
