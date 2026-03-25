import {PlateDimension} from '../components/models'

export const positionFromRowCol = (row: number, col: number, dimension: PlateDimension) => {
  return row * dimension.cols + col
}

export const hrPositionFromPosition = (position: number, dimension: PlateDimension) => {
  const row = Math.floor(position / dimension.cols)
  const col = position - row * dimension.cols + 1
  return `${intToChar(row)}${col}`
}

const intToChar = (int: number) => {
  const code = 'A'.charCodeAt(0)
  return String.fromCharCode(code + int)
}
