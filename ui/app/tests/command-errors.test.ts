/**
 * Tests for the error summary above the command log: the errors must be
 * readable without scrolling through the log.
 */

import { describe, expect, it } from 'vitest'

import type { CommandMessage } from '~/types/management'
import { SHOWN_ERRORS, summarizeCommandErrors } from '~/utils/commandErrors'

const error = (text: string): CommandMessage => ({ level: 'error', text })

describe('summarizeCommandErrors', () => {
  it('shows the error texts, not the info, warnings or the last "Command failed." line', () => {
    const messages: CommandMessage[] = [
      { level: 'info', text: 'Processing file /data/a.csv...' },
      { level: 'warning', text: '64 transfers were not mapped' },
      error('SRC_A -> DST_1: this report was already mapped on 12.09.2026 14:30'),
      error('Command failed.'),
    ]

    expect(summarizeCommandErrors(messages)).toEqual({
      shown: ['SRC_A -> DST_1: this report was already mapped on 12.09.2026 14:30'],
      notShown: 0,
    })
  })

  it('counts the errors it does not show', () => {
    const messages = Array.from({ length: SHOWN_ERRORS + 2 }, (_, index) => error(`Error ${index + 1}`))

    const summary = summarizeCommandErrors(messages)

    expect(summary.shown).toEqual(['Error 1', 'Error 2', 'Error 3', 'Error 4', 'Error 5'])
    expect(summary.notShown).toBe(2)
  })

  it('has nothing to show for a command without errors', () => {
    expect(summarizeCommandErrors([{ level: 'info', text: 'Command completed.' }])).toEqual({
      shown: [],
      notShown: 0,
    })
  })
})
