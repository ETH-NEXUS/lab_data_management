export type NavigationTreeNode = {
  label: string
  icon?: string
  defaultExpanded?: boolean
  children?: NavigationTreeNode[]
  onSelect?: () => void
}
