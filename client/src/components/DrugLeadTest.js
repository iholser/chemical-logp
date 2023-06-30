import React from 'react'

const DrugLeadTest = ({ title, tests, data }) => {
  // { min, max, key }
  const results = tests.map(({ key, min, max }) => {
    let pass = true
    if (min) {
      pass = data[key] >= min
    }
    if (max) {
      pass = pass && data[key] <= max
    }
    return { key, min, max, pass }
  })
  const overall = results.reduce((acc, cur) => acc && cur.pass, true)
  return (
    <div>
      {title} {overall ? 'pass' : 'fail'}
      <div>
        {results.map((r) => (
          <div key={r.key}>
            {r.key}: {r.pass}
          </div>
        ))}
      </div>
    </div>
  )
}

export default DrugLeadTest
