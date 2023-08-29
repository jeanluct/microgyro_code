# waveswim.f90

This repository contains Fortran 90 code for the paper "[Microswimmer trapping in surface waves with shear][1]" by [Francesco Michele Ventrella][2], Nimish Pujara, Guido Boffetta, Massimo Cencini, Jean-Luc Thiffeault and Filippo De Lillo.

## Code and parameter file

The following files part of this repository:

* wavesiwm.f90    - code that integrates the trajectory of a swimmer
* input.dat       - input file

## Citing this work

FM Ventrella, et al. "Microswimmer trapping in surface waves with shear", 	arXiv:2304.14028 [physics.flu-dyn]
DOI: [10.48550/arXiv.2304.14028][3]
BibTeX entry:
```latex
@misc{ventrella2023microswimmer,
      title={Microswimmer trapping in surface waves with shear}, 
      author={Francesco Michele Ventrella and Nimish Pujara and Guido Boffetta and Massimo Cencini and Jean-Luc Thiffeault and Filippo De Lillo},
      year={2023},
      eprint={2304.14028},
      archivePrefix={arXiv},
      primaryClass={physics.flu-dyn}
      abstract{
      Many species of phytoplankton migrate vertically near the surface of the ocean, either in search of light or nutrients. These motile organisms are affected by ocean waves at the surface. We derive a set of wave-averaged equations to describe the motion of spheroidal microswimmers. We include several possible effects, such as gyrotaxis, settling, and wind-driven shear. In addition to the well- known Stokes drift, the microswimmer orbits depend on their orientation in a way that can lead to trapping at a particular depth; this in turn can affect transport of organisms, and may help explain observed phytoplankton layers in the ocean.
      }
}
```

## License

This code is released under the MIT License.  See the file
[LICENSE.txt][4] for copying permission.

[1]: https://arxiv.org/abs/2304.14028 
[2]: francescomichele.ventrella@unito.it 
[3]: https://doi.org/10.48550/arXiv.2304.14028
[4]: https://opensource.org/license/mit/
