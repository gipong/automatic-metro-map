/**
 * Created by gipong on 2017/6/24.
 *
 * provide an easy way for setting parameters in numeric function
 *
 */

(function(window) {
    'use strict';
    function coordsolver() {
        var Coordsolver = {};

        Coordsolver.setVarNumber = function(number) {
            var vars = [];
            var Dmat = [];
            var dvec = [];
            for(var i=0; i<number; i++) {
                vars.push('v'+i);
                Dmat.push(Array.apply(null, Array(number)).map(Number.prototype.valueOf,0));
                dvec.push(0);
            }
            return [vars, Dmat, dvec];
        }

        Coordsolver.addSoftConstraint = function(Dmat, dvec, equation, weight) {
            var regex = /vars\[(\d*)\]/g;
            var coefficient = Array.apply(null, new Array(Dmat[0].length)).map(Number.prototype.valueOf,0);
            var tempLastIndex = 0;
            var constant = 0;

            if(typeof weight == 'undefined') {
                weight = 1;
            }

            var m;
            while((m = regex.exec(equation)) !== null) {
                if (m.index === regex.lastIndex) {
                    regex.lastIndex++;
                }

                coefficient[m[1]] = equation.slice(tempLastIndex, m.index);
                tempLastIndex = regex.lastIndex;
                if(coefficient[m[1]] == '-' || coefficient[m[1]] == '+' || coefficient[m[1]] == '') {
                    coefficient[m[1]] = coefficient[m[1]]+1;
                }
            }

            if(coefficient.length != Dmat[0].length) {
                for(var i=0; i<Dmat[0].length; i++) {
                    if(typeof coefficient[i] == 'undefined') {
                        coefficient[i] = 0;
                    }
                }
            }

            constant = +equation.slice(tempLastIndex);

            coefficient.forEach(function(ele, index) {
                dvec[index] = dvec[index]-weight*ele*constant;
                coefficient.forEach(function(multiele, i) {
                    Dmat[index][i] = Dmat[index][i]+weight*ele*multiele;
                })
            });

            return coefficient;
        }

        Coordsolver.addHardConstraint = function(Amat, bvec, Dmat, equation) {
            var regex = /vars\[(\d*)\]/g;
            var coefficient = Array.apply(null, new Array(Dmat[0].length)).map(Number.prototype.valueOf,0);
            var tempLastIndex = 0;
            var constant = equation.split('>')[1];
            equation = equation.split('>')[0];

            var m;
            while((m = regex.exec(equation)) !== null) {
                if (m.index === regex.lastIndex) {
                    regex.lastIndex++;
                }

                coefficient[m[1]] = equation.slice(tempLastIndex, m.index);
                tempLastIndex = regex.lastIndex;
                if(coefficient[m[1]] == '-' || coefficient[m[1]] == '+' || coefficient[m[1]] == '') {
                    coefficient[m[1]] = coefficient[m[1]]+1;
                }
            }

            if(coefficient.length != Dmat[0].length) {
                for(var i=0; i<Dmat[0].length; i++) {
                    if(typeof coefficient[i] == 'undefined') {
                        coefficient[i] = 0;
                    }
                }
            }

            Amat.push(coefficient);
            bvec.push(constant);

            return [Amat, bvec];
        }

        Coordsolver.printArray = function(array) {
            array.forEach(function(e) {
                console.log(e);
            });
        }

        return Coordsolver;
    }

    if(typeof(Coordsolver) == 'undefined') {
        window.Coordsolver = coordsolver();
    } else {
        console.log('Coordsolver loaded.')
    }
})(window);
