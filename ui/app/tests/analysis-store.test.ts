/**
 * Tests for the analysis store: how the experiment page follows a statistical
 * analysis that runs in the celery-analysis container.
 *
 * The requests to the backend are replaced, so these tests only check what the
 * store does with the answers: which messages it shows, when it gives up, how it
 * comes back after a page reload, and which result list it keeps.
 */

import { createPinia, setActivePinia } from 'pinia'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'

import { useAnalysisStore } from '~/stores/analysis'
import type { LongPollingResponse } from '~/types/management'

const requestApiData = vi.fn()
const requestApiVoid = vi.fn()

vi.mock('~/utils/apiRequests', () => ({
  requestApiData: (...args: unknown[]) => requestApiData(...args),
  requestApiVoid: (...args: unknown[]) => requestApiVoid(...args),
}))

const PAYLOAD = {
  experiment_id: 105,
  label: 'Lum',
  analysis_type: 'single' as const,
  settings: {},
  positive_control: 'P',
  negative_control: 'N',
}
const MINUTE = 60 * 1000

/** The answers of `long_polling`, one per request, the last one repeats. */
const answerWith = (...answers: LongPollingResponse[]) => {
  let index = 0
  requestApiData.mockImplementation((endpoint: string) => {
    if (endpoint.startsWith('long_polling/')) {
      const answer = answers[Math.min(index, answers.length - 1)]
      index += 1
      return Promise.resolve(answer)
    }
    return Promise.resolve({ results: [] })
  })
}

const running = (...texts: string[]): LongPollingResponse => ({
  messages: texts.map((text) => ({ level: 'info', text })),
  next: texts.length,
  status: 'running',
})

// Started, but no worker has taken it yet: no message at all
const waiting: LongPollingResponse = { messages: [], next: 0, status: 'running' }

const completed: LongPollingResponse = {
  messages: [{ level: 'success', text: 'The analysis is done: run.zip' }],
  next: 1,
  status: 'completed',
}

describe('the analysis store', () => {
  beforeEach(() => {
    setActivePinia(createPinia())
    vi.useFakeTimers()
    requestApiData.mockReset()
    requestApiVoid.mockReset().mockResolvedValue(undefined)
    localStorage.clear()
  })

  afterEach(() => {
    vi.useRealTimers()
  })

  it('shows the messages and frees the button when the analysis is done', async () => {
    const store = useAnalysisStore()
    answerWith(running('Running the statistical analysis'), completed)

    const run = store.startAnalysis(PAYLOAD)
    await vi.advanceTimersByTimeAsync(4000)
    await run

    expect(store.messages.map((message) => message.text)).toEqual([
      'Running the statistical analysis',
      'The analysis is done: run.zip',
    ])
    expect(store.status).toBe('completed')
    expect(store.isRunning).toBe(false)
    expect(store.experimentId).toBe(105)
  })

  it('still waits after 10 minutes: the worker may start the analysis until then', async () => {
    const store = useAnalysisStore()
    answerWith(waiting)

    void store.startAnalysis(PAYLOAD)
    await vi.advanceTimersByTimeAsync(10.5 * MINUTE)

    expect(store.isRunning).toBe(true)
    expect(store.messages.some((message) => message.text.startsWith('The analysis has not started yet'))).toBe(true)
  })

  it('gives up after 11 minutes without a message, frees the button and forgets the run', async () => {
    const store = useAnalysisStore()
    answerWith(waiting)

    const run = store.startAnalysis(PAYLOAD)
    await vi.advanceTimersByTimeAsync(11.5 * MINUTE)
    await run

    expect(store.status).toBe('failed')
    expect(store.isRunning).toBe(false)
    expect(store.messages.at(-1)?.text).toContain('did not start within 11 minutes')
    expect(localStorage.getItem('analysis_room_name')).toBeNull()
  })

  it('shows the analysis of this browser again after a page reload', async () => {
    const store = useAnalysisStore()
    localStorage.setItem('analysis_room_name', `105_${Date.now() - MINUTE}`)
    answerWith(running('Running the statistical analysis', 'Step 1 of 3'), completed)

    const resumed = store.resumeAnalysis()
    await vi.advanceTimersByTimeAsync(4000)
    await resumed

    expect(store.experimentId).toBe(105)
    expect(store.messages.map((message) => message.text)).toEqual([
      'Running the statistical analysis',
      'Step 1 of 3',
      'The analysis is done: run.zip',
    ])
    expect(store.status).toBe('completed')
  })

  it('gives up at once on a remembered run that never started', async () => {
    const store = useAnalysisStore()
    localStorage.setItem('analysis_room_name', `105_${Date.now() - 12 * MINUTE}`)
    answerWith(waiting)

    const resumed = store.resumeAnalysis()
    await vi.advanceTimersByTimeAsync(4000)
    await resumed

    expect(store.status).toBe('failed')
    expect(store.isRunning).toBe(false)
  })

  it('keeps the results of the experiment that was opened last, not a late answer', async () => {
    const store = useAnalysisStore()
    let answerFirstExperiment: (value: { results: string[] }) => void = () => {}
    requestApiData.mockImplementation((_endpoint: string, options: { params: { experiment_id: string } }) => {
      if (options.params.experiment_id === '1') {
        return new Promise((resolve) => {
          answerFirstExperiment = resolve
        })
      }
      return Promise.resolve({ results: ['second.zip'] })
    })

    const first = store.fetchResults(1)
    await store.fetchResults(2)
    answerFirstExperiment({ results: ['first.zip'] })
    await first

    expect(store.results).toEqual(['second.zip'])
  })
})
