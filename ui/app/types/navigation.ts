export type NavigationTreeNode = {
  label: string
  icon?: string
  slot?: string
  defaultExpanded?: boolean
  children?: NavigationTreeNode[]
  onSelect?: () => void
  [key: string]: unknown
}
