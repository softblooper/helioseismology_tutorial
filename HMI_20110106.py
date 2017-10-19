import astropy.units as u
from sunpy.net import Fido, attrs as a
from sunpy.database import Database

db = Database()
db.fetch(a.Time("2011-01-06T00:00:00", "2011-01-06T00:01:00"),
         a.Instrument('HMI'), a.vso.Sample(1*u.min))
db.comit()
db
