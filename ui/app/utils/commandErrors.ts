import type { CommandMessage } from '~/types/management'

// The backend ends a failed command with this line. The summary already says
// that the command failed, so it is not repeated there.
const COMMAND_FAILED_LINE = 'Command failed.'

// More errors than this are only counted in the summary; the log shows all of them.
export const SHOWN_ERRORS = 5

/**
 * The errors of a command, for the summary above the log, so that they can be
 * read without scrolling through the log.
 *
 * Accepted data example:
 * - `[{ level: 'info', text: 'Processing file /data/a.csv...' },
 *     { level: 'error', text: 'SRC_A -> DST_1: this report was already mapped on 12.09.2026 14:30, ...' },
 *     { level: 'error', text: 'Command failed.' }]`
 *
 * Returned data example:
 * - `{ shown: ['SRC_A -> DST_1: this report was already mapped on 12.09.2026 14:30, ...'], notShown: 0 }`
 */
export const summarizeCommandErrors = (messages: CommandMessage[]): { shown: string[]; notShown: number } => {
  const errors = messages
    .filter((message) => message.level === 'error')
    .map((message) => message.text)
    .filter((text) => text !== COMMAND_FAILED_LINE)

  return {
    shown: errors.slice(0, SHOWN_ERRORS),
    notShown: Math.max(errors.length - SHOWN_ERRORS, 0),
  }
}
