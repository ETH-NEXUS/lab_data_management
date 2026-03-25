export type PaginatedResponse<T> = {
  results?: T[]
  next?: string | null
}
