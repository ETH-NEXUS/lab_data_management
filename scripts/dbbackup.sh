#!/bin/bash

DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LOG_FILE="${DIR}/backup.log"
BACKUP_DIR="/vol/backups"
FILE_PATTERN="*.backup"

function log {
    echo "$(date +'%Y-%m-%d %H:%M:%S') $@" | tee -a "$LOG_FILE"
}

function db_backup {
    log "Starting DB Backup..."
    cd "${DIR}/.." && docker-compose exec api /bin/bash -c './manage.py dbbackup' 2>&1
    log "DONE DB Backup."
}

function delete_old_backups {
    log "Deleting backups older than 30 days inside Docker..."
    docker-compose exec  api /bin/bash -c "find $BACKUP_DIR -type f -name \"$FILE_PATTERN\" -mtime +30 -exec rm {} +"
    log "Old backups deleted inside Docker."
}

db_backup
delete_old_backups
