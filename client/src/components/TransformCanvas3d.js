import PropTypes from 'prop-types';
import React, { useEffect, useState } from 'react';

import useChemDoodle from '../hooks/useChemDoodle';

const TransformCanvas3d = ({ mol }) => {
  const ChemDoodle = useChemDoodle();
  const [id] = useState(`ChemDoodle-transform-canvas-${Math.floor(Math.random() * 10000)}`);

  useEffect(() => {
    if (ChemDoodle) {
      const transformBallAndStick = new ChemDoodle.TransformCanvas3D(id, 250, 250);
      transformBallAndStick.styles.set3DRepresentation('Ball and Stick');
      transformBallAndStick.styles.backgroundColor = 'white';
      const molecule = ChemDoodle.readMOL(mol, 1);
      transformBallAndStick.loadMolecule(molecule);
    }
  }, [ChemDoodle, id, mol]);

  return (
    <canvas
      className="ChemDoodleWebComponent"
      id={id}
      width="500"
      height="500"
      alt="ChemDoodle Web Component"
      style={{ width: 250, height: 250, backgroundColor: 'white' }}>
      This browser does not support HTML5/Canvas.
    </canvas>
  );
};

TransformCanvas3d.propTypes = {
  mol: PropTypes.string.isRequired,
};

export default TransformCanvas3d;
