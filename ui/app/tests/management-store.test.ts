/**
 * Tests for the management store: how the page follows a running command.
 *
 * The requests to the backend are replaced, so these tests only check what the
 * store does with the answers: which messages it shows, when it stops asking,
 * and when the run button is free again.
 */

import { createPinia, setActivePinia } from 'pinia'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'

import { useManagementStore } from '~/stores/management'
import type { LongPollingResponse } from '~/types/management'

const requestApiData = vi.fn()
const requestApiVoid = vi.fn()

vi.mock('~/utils/apiRequests', () => ({
  requestApiData: (...args: unknown[]) => requestApiData(...args),
  requestApiVoid: (...args: unknown[]) => requestApiVoid(...args),
}))

const FORM_DATA = { command: 'map', machine: 'echo', room_name: 'room_1' }

/** The answers of `long_polling`, one per request, the last one repeats. */
const answerWith = (...answers: LongPollingResponse[]) => {
  let index = 0
  requestApiData.mockImplementation((endpoint: string) => {
    if (endpoint.startsWith('long_polling/')) {
      const answer = answers[Math.min(index, answers.length - 1)]
      index += 1
      return Promise.resolve(answer)
    }
    return Promise.resolve({ directory_content: { type: 'directory', name: '', path: '', children: [] } })
  })
}

const running = (...texts: string[]): LongPollingResponse => ({
  messages: texts.map((text) => ({ level: 'info', text })),
  next: texts.length,
  status: 'running',
})

const finished = (status: 'completed' | 'failed'): LongPollingResponse => ({
  messages: [{ level: 'info', text: `Command ${status}.` }],
  next: 1,
  status,
})

describe('the management store while a command runs', () => {
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

  it('shows the output and frees the button when the command is done', async () => {
    const store = useManagementStore()
    answerWith(running('Processing file a.csv...'), finished('completed'))

    await store.runCommand(FORM_DATA)
    await vi.advanceTimersByTimeAsync(2000)

    expect(store.commandMessages.map((message) => message.text)).toEqual([
      'Executing command: map',
      'Processing file a.csv...',
      'Command completed.',
    ])
    expect(store.commandStatus).toBe('completed')
    expect(store.isRunningCommand).toBe(false)
  })

  it('does not keep the output of an earlier command', async () => {
    const store = useManagementStore()
    answerWith(running('from the first command'))

    await store.runCommand(FORM_DATA)
    await vi.advanceTimersByTimeAsync(1000)
    answerWith(finished('completed'))
    await store.runCommand({ ...FORM_DATA, room_name: 'room_2' })
    await vi.advanceTimersByTimeAsync(3000)

    expect(store.commandMessages.map((message) => message.text)).toEqual([
      'Executing command: map',
      'Command completed.',
    ])
  })

  it('says that a command did not start when nothing is written for a minute', async () => {
    const store = useManagementStore()
    answerWith({ messages: [], next: 0, status: 'running' })

    await store.runCommand(FORM_DATA)
    await vi.advanceTimersByTimeAsync(61_000)

    expect(store.commandStatus).toBe('failed')
    expect(store.isRunningCommand).toBe(false)
    expect(store.commandMessages.at(-1)?.text).toContain('did not start')
  })

  it('asks again when a request for the output fails', async () => {
    const store = useManagementStore()
    requestApiData.mockRejectedValue(new Error('network is down'))

    await store.runCommand(FORM_DATA)
    await vi.advanceTimersByTimeAsync(10_000)

    expect(store.commandStatus).toBe('failed')
    expect(store.isRunningCommand).toBe(false)
    const outputRequests = requestApiData.mock.calls.filter((call) => String(call[0]).startsWith('long_polling/'))
    expect(outputRequests).toHaveLength(5)
  })
})

describe('the management store when the page is opened again', () => {
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

  it('shows the output of the command that is still running', async () => {
    const store = useManagementStore()
    answerWith(running('Processing file a.csv...'), finished('completed'))
    await store.runCommand(FORM_DATA)

    // The page is opened again: a new store with the same browser storage
    setActivePinia(createPinia())
    const newStore = useManagementStore()
    answerWith(running('Processing file a.csv...'), finished('completed'))

    await newStore.resumeCommandOutput()
    await vi.advanceTimersByTimeAsync(2000)

    expect(newStore.commandMessages.map((message) => message.text)).toEqual([
      'Processing file a.csv...',
      'Command completed.',
    ])
    expect(newStore.commandStatus).toBe('completed')
  })

  it('shows the output of a command that ended while the page was closed', async () => {
    const store = useManagementStore()
    answerWith(running('Processing file a.csv...'))
    await store.runCommand(FORM_DATA)

    setActivePinia(createPinia())
    const newStore = useManagementStore()
    answerWith({
      messages: [
        { level: 'error', text: 'The file could not be read.' },
        { level: 'error', text: 'Command failed.' },
      ],
      next: 2,
      status: 'failed',
    })

    await newStore.resumeCommandOutput()

    expect(newStore.commandStatus).toBe('failed')
    expect(newStore.isRunningCommand).toBe(false)
    expect(newStore.commandMessages).toHaveLength(2)
    // The ended command is shown once, the next page load starts empty
    expect(localStorage.getItem('managementRoomName')).toBeNull()
  })

  it('keeps a new command when the answer of the old one comes late', async () => {
    const store = useManagementStore()
    answerWith(running('from the first command'))
    await store.runCommand(FORM_DATA)

    setActivePinia(createPinia())
    const newStore = useManagementStore()
    let answerTheResume: (answer: LongPollingResponse) => void = () => {}
    requestApiData.mockImplementationOnce(() => new Promise((resolve) => (answerTheResume = resolve)))

    const resuming = newStore.resumeCommandOutput()
    answerWith(finished('completed'))
    await newStore.runCommand({ ...FORM_DATA, room_name: 'room_2' })
    answerTheResume(running('from the first command'))
    await resuming
    await vi.advanceTimersByTimeAsync(2000)

    expect(newStore.commandMessages.map((message) => message.text)).toEqual([
      'Executing command: map',
      'Command completed.',
    ])
  })

  it('keeps the command in mind when the first request fails', async () => {
    const store = useManagementStore()
    answerWith(running('Processing file a.csv...'))
    await store.runCommand(FORM_DATA)

    setActivePinia(createPinia())
    const newStore = useManagementStore()
    requestApiData.mockRejectedValue(new Error('network is down'))

    await newStore.resumeCommandOutput()

    expect(localStorage.getItem('managementRoomName')).toBe('room_1')
    expect(newStore.isRunningCommand).toBe(false)
  })
})
