import { Menu, Transition } from '@headlessui/react'
import React, { Fragment, useState } from 'react'
import Sketcher from './components/Sketcher'
import TransformCanvas3d from './components/TransformCanvas3d'

import { withChemDoodle } from './hooks/useChemDoodle'

const examples = [
  {
    name: 'Fluoxetine',
    smiles: 'CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F',
  },
  {
    name: 'Benzene',
    smiles: 'c1ccccc1',
  },
  {
    name: 'Anthracene',
    smiles: 'C1=CC=C2C=C3C=CC=CC3=CC2=C1',
  },
  {
    name: 'Aflatoxin B2',
    smiles: 'COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5CCO[C@@H]5OC4=C1',
  },
  {
    name: 'Aspirin',
    smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O',
  },
]

const url = 'https://chem-logp-2xip2mubea-uw.a.run.app'
// const url = 'http://localhost:8000'

function classNames(...classes) {
  return classes.filter(Boolean).join(' ')
}

const Chem = () => {
  const [loading, setLoading] = useState(false)
  const [results, setResults] = useState()

  const getResultsMol = async (mol) => {
    setLoading(true)
    const response = await fetch(`${url}/mol`, {
      method: 'POST',
      body: mol,
    }).catch((err) => {
      console.log(err)
      return { json: () => ({ error: 'Failed to predict chemical' }) }
    })
    setLoading(false)
    if (response.status === 500) {
      console.log(response)
      setResults({ error: 'Failed to predict chemical' })
    } else {
      setResults(await response.json())
    }
  }

  const getResultsSmiles = async (smiles) => {
    setLoading(true)
    const response = await fetch(`${url}/smiles`, {
      method: 'POST',
      body: smiles,
    })
    setLoading(false)
    setResults(await response.json())
  }

  let percentOctanol = 100
  if (results?.logP) {
    percentOctanol = 100 * (10 ** results.logP / (10 ** results.logP + 1))
  }

  return (
    <div className="m-2">
      {results?.error && <div>An Error Occured</div>}
      {loading && <div>loading...</div>}
      {results?.mol ? (
        <div>
          <button
            onClick={() => setResults(null)}
            className="p-2 bg-blue-200 hover:bg-blue-500 rounded"
          >
            &larr;clear
          </button>
          <div className="grid grid-rows-2 grid-cols-2 gap-2 bg-white">
            <div className="rounded overflow-hidden shadow-lg m-1 p-2 grid justify-items-center">
              <TransformCanvas3d mol={results.mol} />
            </div>
            <div className="rounded overflow-hidden shadow-lg m-1 p-2">
              <div>
                <span className="text-2xl">logP: </span>
                <span className="text-3xl font-bold">{results?.logP}</span>
              </div>
              <div>{percentOctanol.toFixed(2)}% in octanol</div>
              <div className="h-1/2">
                <div className="h-full bg-blue-400 w-5 mb-6">
                  <div
                    className="bg-orange-200 w-5"
                    style={{ height: `${percentOctanol.toFixed(2)}%` }}
                  ></div>
                </div>
                <div>{(100 - percentOctanol).toFixed(2)}% in water</div>
              </div>
            </div>
            <div>
              <h2 className="text-2xl">Interpretation</h2>
              {results.logP <= 5 && results.logP >= -0.4 ? (
                <div>likely drug candidate (logP between -0.4 and 5)</div>
              ) : (
                <div>
                  not bioavailable, unlikely to have potential as a drug
                </div>
              )}
              <div>
                <h2 className="text-xl">Environmental Fate</h2>
                {results.logP < 0 && (
                  <div>Potential for groundwater contamination</div>
                )}
                {results.logP > 2 && (
                  <div>Likey to cause soil contamination.</div>
                )}
                {results.logP > 4 && <div>Likey to bio accumulate.</div>}
              </div>
            </div>
            <div>
              Performance
              <div>
                <h2>{results.logP}</h2>
              </div>
              <div>Wildman Crippen: {results.wildman_crippen_logP}</div>
              <div>Experimental: {results.experimental_logp}</div>
              <div>XLogP3: {results.XLogP3}</div>
            </div>
          </div>
        </div>
      ) : (
        <div className="grid justify-items-center bg-white">
          <Sketcher onSubmit={getResultsMol} submitLabel="Predict&rarr;" />

          <Menu as="div" className="relative inline-block text-left">
            <div>
              <Menu.Button className="inline-flex justify-center w-full rounded-md border border-gray-300 shadow-sm px-4 py-2 bg-white text-sm font-medium text-gray-700 hover:bg-gray-50 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-offset-gray-100 focus:ring-indigo-500">
                Example Chemicals
              </Menu.Button>
            </div>

            <Transition
              as={Fragment}
              enter="transition ease-out duration-100"
              enterFrom="transform opacity-0 scale-95"
              enterTo="transform opacity-100 scale-100"
              leave="transition ease-in duration-75"
              leaveFrom="transform opacity-100 scale-100"
              leaveTo="transform opacity-0 scale-95"
            >
              <Menu.Items className="origin-top-right absolute right-0 mt-2 w-56 rounded-md shadow-lg bg-white ring-1 ring-black ring-opacity-5 focus:outline-none">
                <div className="py-1">
                  {examples.map((ex) => (
                    <Menu.Item key={ex.name}>
                      {({ active }) => (
                        <button
                          onClick={() => getResultsSmiles(ex.smiles)}
                          className={classNames(
                            active
                              ? 'bg-gray-100 text-gray-900'
                              : 'text-gray-700',
                            'block px-4 py-2 text-sm w-full',
                          )}
                        >
                          {ex.name}
                        </button>
                      )}
                    </Menu.Item>
                  ))}
                </div>
              </Menu.Items>
            </Transition>
          </Menu>
        </div>
      )}
    </div>
  )
}

export default withChemDoodle(Chem)
