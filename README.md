# Unsaturated_Flow_2D_Steady-state

![Preview](https://numericalenvironmental.files.wordpress.com/2018/05/sat_script.png?w=1632)

This script, written in Julia (v. 1.0 and above), solves the steady-state form of the Richards equation for variably-saturated flow in porous media in a two-dimensional vertical cross section. Model output can be used, for example, to understand the distribution of soil moisture or the formation of perched aquifers in heterogeneous soils under different types of boundary conditions. The model iteratively solves the non-linear flow problem by solving a matrix of flow-balance equations for recent-estimate values of matrix potential, then updating hydraulic conductivity values as a function of matrix potential. More information is available in my blog post, https://numericalenvironmental.wordpress.com/2018/05/26/a-steady-state-variably-saturated-flow-model-in-vertical-cross-section-a-finite-difference-approach-using-julia/.

The following space- or tab-delimited text input files are all required:

* Materials.txt - geologic material properties table (name, horizontal and vertical saturated hydraulic conductivities, Van Genuchten alpha and N parameters, residual saturation)
* Domain.txt - spatial extent and discretization of model; specify default geologic material 
* Blocks.txt - list of rectangular heterogeneities (bottom left and top right coordinates, corresponding geologic material name)
* Cells.txt - an alternative means of specifying heterogeneities, assigning properties and source fluxes on a cell-by-cell basis
* Sources.txt - list of line boundary conditions (type - whether flux or specified pressure, orientation of line - x- or z-directions, position of line with reference to the perpendicular axis, starting and ending points along its parallel axis, and value - whether inflowing velocity or fixed pressure)
* Solver.txt - initial uniform value for pressure head as a guess (to inform hydraulic conductivity estimates), iteration result weighting factor (see link to blog post, above, for details), and number of iterations

Model output, which includes flux vectors and flow balances, is written to a comma-delimited text file.

Note that the example problem described by the input files is set to read hydraulic properties from the Cells.txt file, as opposed to the Blocks.txt and Sources.txt files. The latter are provided here for reference only to illustrate formatting.

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

