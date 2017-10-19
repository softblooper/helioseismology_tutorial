import astropy.units as u
from sunpy.net import Fido, attrs as a
from sunpy.database import Database

db = Database()
db.fetch(a.Time("2011-01-05T23:55:00", "2011-01-05T23:56:00"),
         a.Instrument('MDI'), a.vso.Sample(1*u.min))
db.comit()
db
