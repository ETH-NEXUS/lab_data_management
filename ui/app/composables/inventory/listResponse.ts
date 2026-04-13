import type { PaginatedResponse } from '~/types/api'
import { toRelativeApiEndpoint } from '~/utils/apiEndpoint'

export type ListResponse<T> = PaginatedResponse<T> | T[]

export type ParsedListPage<T> = {
  items: T[]
  nextEndpoint: string | null
}

/**
 * Normalizes list endpoint payloads to one common shape.
 *
 * Accepted response examples:
 * - paginated: `{ results: [{ id: 1 }], next: '/api/inventory/stocks/?page=2' }`
 * - plain list: `[{ id: 1 }, { id: 2 }]`
 *
 * Returned data examples:
 * - `{ items: [{ id: 1 }], nextEndpoint: 'inventory/stocks/?page=2' }`
 * - `{ items: [{ id: 1 }, { id: 2 }], nextEndpoint: null }`
 */
export const parseListPage = <T>(payload: ListResponse<T>): ParsedListPage<T> => {
  if (Array.isArray(payload)) {
    return {
      items: payload,
      nextEndpoint: null,
    }
  }

  return {
    items: payload.results ?? [],
    nextEndpoint: toRelativeApiEndpoint(payload.next ?? null),
  }
}
