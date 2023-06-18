# fermi2

Code supporting the article [Gamma-Ray Bubbles on the Sea of Wavelets: Quantifying Existence and Morphology of the Fermi Bubbles](https://storage.googleapis.com/fermi-tipsh/title.html).

This is kind of a mess, and some day I will try to clean it up. The notebooks do run,
but you will need the data, which I have to clean up and upload someplace. Ditto the
various python scripts in the top-level folder, they need the processed Fermi data, which
is too large to put in github. The core algorithm code is in the `fermi` folder, and should
be solid. There's some nice tricks in there to make it go fast. The one constraint is that
it's oriented for all-sky analysis, so the data is reflected at the boundaries assuming
it is projected on the entire celestial sphere.

# LICENSE

Copyright © 2023 David D. Dixon

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
