import React, { useEffect, useState } from 'react'

const createElement = async (type, props) => {
  try {
    const el = document.createElement(type)
    Object.assign(el, props)
    document.body.appendChild(el)
    await new Promise((r) => setTimeout(r, 100))
  } catch (e) {
    console.log(e)
  }
}

const getChemDoodle = async () => {
  return new Promise((resolve) => {
    const i = setInterval(() => {
      if (window.ChemDoodle) {
        clearInterval(i)
        resolve(window.ChemDoodle)
      }
    }, 50)
  })
}

const initChemDoodle = async () => {
  // await new Promise((r) => setTimeout(r, 50))
  // if (!window.ChemDoodle && !window.initChemDoodle) {
  //   window.initChemDoodle = true
  //   await createElement('link', {
  //     href: '/lib/ChemDoodle-ui.css',
  //     rel: 'stylesheet',
  //   })
  //   await createElement('script', { async: true, src: '/lib/ChemDoodleWeb.js' })
  //   await createElement('script', {
  //     async: true,
  //     src: '/lib/ChemDoodleWeb-uis.js',
  //   })
  //   await createElement('script', {
  //     innerText: 'window.ChemDoodle = ChemDoodle;',
  //   })
  // } else {
  //   await new Promise((r) => setTimeout(r, 10))
  // }

  const chemdoodle = await getChemDoodle()
  return chemdoodle
}

export const withChemDoodle = (Comp) => (props) => {
  const [ChemDoodle, setChemDoodle] = useState(window.ChemDoodle)
  useEffect(() => {
    if (!ChemDoodle) {
      console.log('no chemdoodle')
      initChemDoodle().then((c) => setChemDoodle(c))
    }
  }, [ChemDoodle, setChemDoodle])

  if (!ChemDoodle) return <div>ChemDoodle loading</div>
  return <Comp ChemDoodle={ChemDoodle} {...props} />
}

const useChemDoodle = () => window.ChemDoodle

export default useChemDoodle
