/* eslint-disable max-len */
import gql from 'graphql-tag';

const inventoryReferenceFields = `{
  _id
  name
  campusCode
  tenantCode
}`;

const locationReferenceFields = `{
  buildingId
  floorId
  roomId
  sublocationId
  buildingName
  floorName
  roomNumber
  sublocationName
}`;

export const familyFields = gql`
  fragment familyFields on Family {
    _id
    identifiers {
      cas
    }
    components {
      cas
    }
    boilingPoint
    density
    description
    flags {
      id
      date
      type
      problem
    }
    flashPoint
    formula
    form
    formNormalize
    meltingPoint
    molecularWeight
    name
    scope
    smiles
    synonyms
    firstAid {
      general
      inhalation
      skin
      eye
      ingestion
      fire
      explosion
      exposure
    }
    storage {
      protection
      requirements
      remarks
    }
    healthSymptoms {
      general
      inhalation
      skin
      eye
      ingestion
      fire
      explosion
      exposure
    }
    nfpa {
      flammability
      health
      reactivity
      specialHazard
    }
    ghs {
      hazardCode
      hazardClass
      category
      hazardStatement
      signalWord
    }
    pictograms
    sds {
      particleSize
    }
  }
`;

export const CONTAINER_FRAGMENT = gql`fragment containerFields on Container {
  _id
  audits {
    type,
    rule,
    message,
    actionMessage,
    field,
    value,
    fieldType,
    ignoreBy {
      userId
      firstName
      lastName
    }
  }
  barcode,
  customBarcode
  size
  amount
  units
  type
  orderId
  orderItemId
  openedDate
  receivedDate
  lastTestedDate
  comments
  isPrivate
  tags
  expirationDate
  physicalState
  concentration
  concentrationUnits
  solvent
  checkout
  transfer {
    name
    targetInventory {
      _id
      name
      campusCode
      tenantCode
    }
    initiatedBy {
      userId
      firstName
      lastName
      tenantCode
      campusCode
    }
  }
  family {
    _id
    density
    name
  }
  substance {
    _id
  },
  inventory ${inventoryReferenceFields}
  location ${locationReferenceFields}
}`;

export const FAMILY_SEARCH_RESULT_FIELDS = gql`
  fragment familySearchResultFields on FamilySearchResult {
    _id
    identifiers {
      cas
    }
    name
    components
    form
    formNormalize
    formula
    ghs {
      hazardCode
    }
    pictograms
    tags
    inventoryNames
    totalContainers
    bands
    audits
    highlight
  }
`;

export const FAMILY_COMMON_FIELDS = `
  ...familyFields
  substances {
    _id
    name
    vendor
    products {
      productNumber
    }
  }
  inventoryAttachments(inventoryId: $inventoryId) {
    _id
    attachmentType
    attachment {
      id
      filename
      token
    }
  }
  customName(inventoryId: $inventoryId) {
    name
  }
  containers(campusCode: $campusCode, inventoryIds: [$inventoryId], page: $page, displayCount: $displayCount,
    barcode: $barcode, sublocationId: $sublocationId, roomId: $roomId, sortCriteria: $sortCriteria) {
    ...containerFields
  }
  attachments {
    id
    filename
    token
  }
  bands {
    name
    class
    requireSop
    campus
    useRules
  }
  containerCount(inventoryIds: [$inventoryId])
  sops(inventoryId: $inventoryId) {
    _id
    document {
      title
      description
      label
      acknowledgement {
        showAcknowledgement
        requiredSignatures
        signatures
      }
    }
  }
`;
