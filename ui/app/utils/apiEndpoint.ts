/**
 * Converts API pagination links into the relative endpoint format expected by useAPI().
 *
 * Accepted input examples:
 * - `null`
 * - `/api/projects/?page=2`
 * - `https://host/api/experiments/?page=3`
 *
 * Returned value examples:
 * - `null`
 * - `projects/?page=2`
 * - `experiments/?page=3`
 */
export const toRelativeApiEndpoint = (input: string | null | undefined): string | null => {
  if (!input) return null

  let endpoint = input
  if (endpoint.startsWith('http://') || endpoint.startsWith('https://')) {
    try {
      const parsed = new URL(endpoint)
      endpoint = `${parsed.pathname}${parsed.search}`
    } catch {
      return null
    }
  }

  return endpoint.replace(/^\/api\//, '').replace(/^\//, '')
}
