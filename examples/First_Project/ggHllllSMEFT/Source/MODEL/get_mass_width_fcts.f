      DOUBLE PRECISION FUNCTION GET_MASS_FROM_ID(ID)
      IMPLICIT NONE
      INTEGER ID
      INCLUDE 'coupl.inc'

      IF (ID.EQ.23) THEN
        GET_MASS_FROM_ID=ABS(MDL_MZ)
      ELSE IF (ID.EQ.24.OR.ID.EQ.-24) THEN
        GET_MASS_FROM_ID=ABS(MDL_MW)
      ELSE IF (ID.EQ.9000002.OR.ID.EQ.-9000002) THEN
        GET_MASS_FROM_ID=ABS(MDL_MZ)
      ELSE IF (ID.EQ.9000003.OR.ID.EQ.-9000003) THEN
        GET_MASS_FROM_ID=ABS(MDL_MW)
      ELSE IF (ID.EQ.9000004.OR.ID.EQ.-9000004) THEN
        GET_MASS_FROM_ID=ABS(MDL_MW)
      ELSE IF (ID.EQ.6.OR.ID.EQ.-6) THEN
        GET_MASS_FROM_ID=ABS(MDL_MT)
      ELSE IF (ID.EQ.25) THEN
        GET_MASS_FROM_ID=ABS(MDL_MH)
      ELSE
        GET_MASS_FROM_ID=0D0
      ENDIF
      RETURN
      END


      DOUBLE PRECISION FUNCTION GET_WIDTH_FROM_ID(ID)
      IMPLICIT NONE
      INTEGER ID
      INCLUDE 'coupl.inc'

      IF (ID.EQ.23) THEN
        GET_WIDTH_FROM_ID=ABS(MDL_WZ)
      ELSE IF (ID.EQ.24.OR.ID.EQ.-24) THEN
        GET_WIDTH_FROM_ID=ABS(MDL_WW)
      ELSE IF (ID.EQ.9000002.OR.ID.EQ.-9000002) THEN
        GET_WIDTH_FROM_ID=ABS(MDL_WZ)
      ELSE IF (ID.EQ.9000003.OR.ID.EQ.-9000003) THEN
        GET_WIDTH_FROM_ID=ABS(MDL_WW)
      ELSE IF (ID.EQ.9000004.OR.ID.EQ.-9000004) THEN
        GET_WIDTH_FROM_ID=ABS(MDL_WW)
      ELSE IF (ID.EQ.6.OR.ID.EQ.-6) THEN
        GET_WIDTH_FROM_ID=ABS(MDL_WT)
      ELSE IF (ID.EQ.25) THEN
        GET_WIDTH_FROM_ID=ABS(MDL_WH)
      ELSE
        GET_WIDTH_FROM_ID=0D0
      ENDIF
      RETURN
      END


