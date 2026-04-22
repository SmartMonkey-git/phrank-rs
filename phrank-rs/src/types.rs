pub(crate) type PatientPhenotypePair<'a> = (&'a String, &'a Vec<String>);
pub(crate) type PatientPhenotypePairRef<'a> =
    &'a (PatientPhenotypePair<'a>, PatientPhenotypePair<'a>);
