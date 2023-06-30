import PropTypes from 'prop-types'
import React, { useEffect, useState } from 'react'

import { withChemDoodle } from '../hooks/useChemDoodle'

const Sketcher = ({ ChemDoodle, mol, onChange, onSubmit, submitLabel }) => {
  const [id] = useState(
    `ChemDoodle-sketcher-${Math.floor(Math.random() * 10000)}`,
  )
  const [sketcher, setSketcher] = useState()

  useEffect(() => {
    const sketch = new ChemDoodle.SketcherCanvas(id, 500, 300, {
      useServices: false,
      oneMolecule: true,
    })
    if (mol) {
      const molecule = ChemDoodle.readMOL(mol, 30)
      sketch.loadMolecule(molecule)
    }
    sketch.checkScale()
    sketch.repaint()
    setSketcher(sketch)
  }, [ChemDoodle, id, mol, setSketcher])

  const getMol = () => onSubmit(ChemDoodle.writeMOL(sketcher.getMolecule()))
  return (
    <div className="grid justify-items-center bg-white shadow-lg p-4">
      <canvas
        className="ChemDoodleWebComponent"
        id={id}
        width="500"
        height="300"
        alt="ChemDoodle Web Component"
        style={{ width: 500, height: 30, backgroundColor: 'white' }}
      >
        This browser does not support HTML5/Canvas.
      </canvas>
      <div className="flex w-full justify-end">
        <button
          onClick={getMol}
          className="p-2 bg-blue-200 hover:bg-blue-500 rounded"
        >
          {submitLabel}
        </button>
      </div>
    </div>
  )
}

Sketcher.propTypes = {
  mol: PropTypes.string,
  onChange: PropTypes.func,
  onSubmit: PropTypes.func,
  submitLabel: PropTypes.string,
}

Sketcher.defaultProps = {
  mol: null,
  onChange: () => {},
  onSubmit: console.log,
  submitLabel: 'Submit',
}

export default withChemDoodle(Sketcher)
