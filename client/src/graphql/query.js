/* eslint-disable max-len */
import gql from 'graphql-tag';

import { CONTAINER_FRAGMENT, familyFields, FAMILY_COMMON_FIELDS, FAMILY_SEARCH_RESULT_FIELDS } from './fragment';

export const INVENTORY_FAMILY_DETAIL = gql`
  query LoadInventoryDetail($id: ID!, $campusCode: String, $inventoryId: String, $page: Int, $displayCount: Int, $barcode: ID, $sublocationId: ID, $roomId: ID, $sortCriteria: String) {
    family(id: $id, inventoryId: $inventoryId) {
      ${FAMILY_COMMON_FIELDS}
    }
  }
  ${familyFields}
  ${CONTAINER_FRAGMENT}
`;

export const SEARCH_FAMILIES = gql`
  query SearchFamily($campusCode: String!, $keyword: String) {
    families(campusCode: $campusCode, keyword: $keyword) {
      ...familySearchResultFields
    }
  }
  ${FAMILY_SEARCH_RESULT_FIELDS}
`;
