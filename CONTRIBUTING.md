# Contributing to GaiaVario

Thank you for considering a contribution. GaiaVario is developed at the
Deforming Gaia Lab, Department of Geology, University of Patras, within the
TRANSFORM² project.

## Before you start

- Open an issue describing the bug or the change. Changes that alter the
  solution (models, weights, rejection) need a before/after comparison on a
  recorded stream and the synthetic-truth tests kept green.
- Never include caster credentials, receiver addresses, user names or station
  configuration files in an issue, a commit or a log excerpt. Scrub
  `config/ntrip_profiles.toml` and session logs before attaching them.

## Pull requests

- One topic per pull request, with a short description of what changed and
  how it was tested (`julia --project=. test/runtests.jl`; for
  changes that alter solutions also a replay comparison on a capture).
- Keep the style of the surrounding code; update `README.md` and
  `CHANGELOG.md` when behaviour changes.

## Licence of contributions and sign-off

GaiaVario is licensed under the GNU Affero General Public License v3.0 only.
By submitting a contribution you certify the
[Developer Certificate of Origin 1.1](https://developercertificate.org/) and
agree that your contribution is licensed under AGPL-3.0-only, and that the
project maintainer may also make it available, together with the rest of the
project, under a later version of the AGPL. Add a sign-off line to each
commit:

```
Signed-off-by: Your Name <your@email>
```

(`git commit -s` adds it.)

## Contact

Panagiotis Elias — <pelias@upatras.gr>
