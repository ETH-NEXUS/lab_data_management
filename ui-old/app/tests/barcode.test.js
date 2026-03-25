import {generateBarcodes} from '../src/components/helpers.ts'
import {describe, expect, it} from 'vitest'

describe('generateBarcodes', () => {
  it('generates barcodes with specified prefix, number of plates, and sides', () => {
    const prefix = 'ABC'
    const numberOfPlates = 3
    const sides = ['North', 'East']
    const expectedBarcodes = [
      {NorthBarcode: 'ABC_1$ABC_1', SouthBarcode: '', EastBarcode: 'ABC_1$ABC_1', WestBarcode: ''},
      {NorthBarcode: 'ABC_2$ABC_2', SouthBarcode: '', EastBarcode: 'ABC_2$ABC_2', WestBarcode: ''},
      {NorthBarcode: 'ABC_3$ABC_3', SouthBarcode: '', EastBarcode: 'ABC_3$ABC_3', WestBarcode: ''},
    ]

    const barcodes = generateBarcodes(prefix, numberOfPlates, sides)

    expect(barcodes).toEqual(expectedBarcodes)
  })

  it('returns empty array if number of plates is 0', () => {
    const prefix = 'ABC'
    const numberOfPlates = 0
    const sides = ['North', 'East']

    const barcodes = generateBarcodes(prefix, numberOfPlates, sides)

    expect(barcodes).toEqual([])
  })
})
