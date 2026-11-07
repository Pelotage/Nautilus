//! Type-safe physical quantities. Nothing else may define a physical quantity.

/// Days per year (Julian convention, 365.25 days).
pub const DAYS_PER_YEAR: f64 = 365.25;

/// Metres per nautical mile (exact by definition).
const METRES_PER_NAUTICAL_MILE: f64 = 1852.0;

/// Seconds per hour.
const SECONDS_PER_HOUR: f64 = 3600.0;

/// Time in days.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Days(pub f64);

impl Days {
    /// Returns the underlying value in days.
    pub fn value(self) -> f64 {
        self.0
    }

    /// Converts a duration in hours into days.
    pub fn from_hours(h: f64) -> Days {
        Days(h / 24.0)
    }

    /// Converts this duration into years, using [`DAYS_PER_YEAR`].
    pub fn years(self) -> f64 {
        self.0 / DAYS_PER_YEAR
    }
}

/// Temperature in degrees Celsius.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Celsius(pub f64);

impl Celsius {
    /// Returns the underlying value in degrees Celsius.
    pub fn value(self) -> f64 {
        self.0
    }
}

/// Salinity in practical salinity units (PSU, dimensionless).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Psu(pub f64);

impl Psu {
    /// Returns the underlying value in practical salinity units.
    pub fn value(self) -> f64 {
        self.0
    }
}

/// Speed in knots (nautical miles per hour).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Knots(pub f64);

impl Knots {
    /// Returns the underlying value in knots.
    pub fn value(self) -> f64 {
        self.0
    }

    /// Converts this speed into metres per second (1 knot = 1852/3600 m/s).
    pub fn to_metres_per_second(self) -> f64 {
        self.0 * METRES_PER_NAUTICAL_MILE / SECONDS_PER_HOUR
    }
}

/// Length in metres.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Metres(pub f64);

impl Metres {
    /// Returns the underlying value in metres.
    pub fn value(self) -> f64 {
        self.0
    }
}

/// Length in kilometres.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Kilometres(pub f64);

impl Kilometres {
    /// Returns the underlying value in kilometres.
    pub fn value(self) -> f64 {
        self.0
    }
}

/// Mass concentration in milligrams per cubic metre (used for chlorophyll-a).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct MilligramsPerCubicMetre(pub f64);

impl MilligramsPerCubicMetre {
    /// Returns the underlying value in mg/m³.
    pub fn value(self) -> f64 {
        self.0
    }
}

/// Length in micrometres (1e-6 m).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Micrometres(pub f64);

impl Micrometres {
    /// Returns the underlying value in micrometres.
    pub fn value(self) -> f64 {
        self.0
    }

    /// Converts this length into metres.
    pub fn to_metres(self) -> f64 {
        self.0 * 1e-6
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    /// Unit conversions: knots → m/s, hours → days, days → years, µm → m.
    #[test]
    fn conversions() {
        assert!((Knots(1.0).to_metres_per_second() - 0.514_444_444_444).abs() < 1e-9);
        assert_eq!(Days::from_hours(12.0).value(), 0.5);
        assert!((Days(365.25).years() - 1.0).abs() < 1e-12);
        assert!((Micrometres(150.0).to_metres() - 1.5e-4).abs() < 1e-15);
        assert_eq!(Celsius(25.0).value(), 25.0);
    }
    /// Chlorophyll-a concentration unit holds its value (mg/m³).
    #[test]
    fn chlorophyll_unit() {
        assert_eq!(MilligramsPerCubicMetre(2.5).value(), 2.5);
    }
}
