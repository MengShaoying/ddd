function min(v){
	for (var i = 1; i < v.length; i++){
		if (v[i] < v[0]) v[0] = v[i];
	}
	return v[0];
}

function max(v){
	for (var i = 1; i < v.length; i++){
		if (v[i] > v[0]) v[0] = v[i];
	}
	return v[0];
}

function draw(x,y){
	var dy = max(y) - min(y);
	var dx = max(x) - min(x);
	var k = 300/dy;
	var h = '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" height="'+k*dy+'" width="'+k*dx+'">';
	for (var i = 0; i < x.length; i++){
		h += '<circle cx="'+k*(x[i]-min(x))+'" cy="'+k*(dy-y[i]+min(y))+'" r="2" fill="black"/>';
	}
	h += '</svg>';
	var p = document.createElement('div');
	p.innerHTML = h;
	document.body.appendChild(p);
}

/*
 * 页面上输出变量
 */
function var_dump(param){
	var p = document.createElement('p');
	switch (typeof(param)){
	case 'string':
		p.innerText = param;
		break;
	case 'number':
		p.innerText = param;
		break;
	case 'undefined':
		p.innerText = 'undefined';
		break;
	default:
		p.innerText = JSON.stringify(param);
	}
	document.body.appendChild(p);
}

/*
 * 高斯函数的Javascript实现
 */
function gaussianFunction(x, mu, sigma){
	var r2 = 0;
	if ('number' == typeof(x) && 'number' == typeof(mu)){
		r2 = (x - mu) * (x - mu);
		return Math.exp((0-r2)/(2*sigma*sigma))/Math.sqrt(2*Math.PI*sigma);
	} else {
		for (var i = 0; i < x.length; i++) r2 += (x[i] - mu[i]) * (x[i] - mu[i]);
		return Math.exp((0-r2)/(2*sigma*sigma))/Math.pow(Math.sqrt(2*Math.PI*sigma),x.length);
	}
}
//var N = 3;
//var x = [3,3,4];
//var u = [5,5,5];
//var sigma = 2;
//var y = gaussianFunction(x, u, sigma, N);

/*
 * 查看矩阵
 */
function viewMatrix(matrix){
	var html = '';
	for (var i = 0; i < matrix.length; i++){
		for (var j = 0; j < matrix[i].length; j++){
			html += matrix[i][j] + '&nbsp;';
		}
		html += '<br/>';
	}
	var p = document.createElement('p');
	p.innerHTML = html;
	document.body.appendChild(p);
}

/*
 * 矩阵乘法
 */
function matrixMultiplication(a,b){
	var result = Array();
	for (var i = 0; i < a.length; i++){
		result.push(Array());
		for (var j = 0; j < b[0].length; j++){
			result[i].push(0);
			for (var k = 0; k < b.length; k++)
				result[i][j] += a[i][k] * b[k][j];
		}
	}
	return result;
}
//var a = [[1,2,3],[4,5,6]];
//var b = [[1,4],[2,5],[3,6]];
//viewMatrix(a);
//viewMatrix(b);
//viewMatrix(matrixMultiplication(a,b));

/*
 * 矩阵加法
 */
function addMatrix(a, b){
	if ('number' == typeof(a)){
		for (var i = 0; i < b.length; i++)
			for (var j = 0; j < b[i].length; j++)
				b[i][j] += a;
		return b;
	}
	for (var i = 0; i < b.length; i++)
		for (var j = 0; j < b[i].length; j++)
			b[i][j] += a[i][j];
	return b;
}

function decMatrix(a, b){
	for (var i = 0; i < b.length; i++)
		for (var j = 0; j < b[i].length; j++)
			b[i][j] -= a[i][j];
	return b;
}

//var a = [[1,2,3],[4,5,6],[7,8,9]];
//var b = [[1,1,1],[1,1,1],[1,1,1]];
//var c = addMatrix(a,b);
//viewMatrix(c);

/*
 * 求已知矩阵的转置矩阵
 */
function matrixT(matrix){
	var r, t;
	t = Array();
	for (var j = 0; j < matrix[0].length; j++){
		r = Array();
		for (var i = 0; i < matrix.length; i++){
			r.push(matrix[i][j]);
		}
		t.push(r);
	}
	return t;
}
//var a = [[1,2,2],[4,5,6]];
//viewMatrix(a);
//viewMatrix(matrixT(a));

/*
 * 计算矩阵的F-范数
 */
function F_Norm(matrix){
	var s = 0;
	for (var i = 0; i < matrix.length; i++)
		for (var j = 0; j < matrix[i][j]; j++)
			s += matrix[i][j] * matrix[i][j];
	return Math.sqrt(s);
}

function F_Norm2(matrix){
	var s = 0;
	for (var i = 0; i < matrix.length; i++)
		for (var j = 0; j < matrix[i][j]; j++)
			s += matrix[i][j] * matrix[i][j];
	return s;
}

/* 求和 */
function E(vector,n){
	var sum = 0;
	for (var i = 0; i < n; i++)
		sum += vector[i];
	return sum;
}
//var a = Array(5,4,2,7,8);
//var_dump(a);
//var_dump(E(a,5));

/*
 * Parameter: matrix must be a invertible matrix.
 * matrix.length = matrix[i].length
 * i=0,1,2,3,....matrix.length.
 * Parameter: matrix looks like:
 * [ [1,2,3,4],
 *   [2,4,3,1],
 *   [1,2,8,5] ]
 */
function getInverseMatrix(matrix){
	for (var i = 0; i < matrix.length; i++){
		for (var j = 0; j < matrix.length; j++)
			matrix[i].push(0);
		matrix[i][i + matrix.length] = 1;
	}
	var diagElement, multi;
	for (var i = 0; i < matrix.length; i++){
		diagElement = matrix[i][i];
		for (var j = 0; j < matrix[0].length; j++)
			matrix[i][j] /= diagElement;
		for (var k = 0; k < matrix.length; k++){
			if (k == i) continue;
			multi = matrix[k][i];
			for (var n = 0; n < matrix[0].length; n++)
				matrix[k][n] -= matrix[i][n] * multi;
		}
	}
	for (var i = 0; i < matrix.length; i++)
		for (var j = 0; j < matrix.length; j++)
			matrix[i].shift();
	return matrix;
}

function diag(vector){
	var matrix = Array();
	for (var i = 0; i < vector.length; i++){
		matrix.push(Array());
		for (var j = 0; j < vector.length; j++)
			matrix[i].push(0);
	}
	for (var i = 0; i < vector.length; i++){
		matrix[i][i] = vector[i];
	}
	return matrix;
}

/* RVM核函数 */
function PHI(x){
	var p = Array();var s = Array();
	for (var i = 0; i < x.length; i++){
		s = Array();s.push(1);
		for (j = 0; j < x.length; j++)
			s.push(gaussianFunction(x[i], x[j], 1));
		p.push(s);
	}
	return p;
}

function RVM(){
	this.x = Array();
	this.t = Array();
	this._alpha = Array();
	this._sigma = 0;
	this._PHI   = Array();
	this._A     = Array();
	this._SIGMA = Array();
	this._u     = Array();
	this._gama  = Array();
	this.setTrainXT = function(x,t){
		this.x  = x;
		this.t  = t;
	};
	this.initAlphaSigma = function(){
		this._alpha = Array();
		for (var i = 0; i < this._N + 1; i++)
			this._alpha.push(1 + Math.random());
		this._sigma = 1 + Math.random();
	};
	this.train = function(){
		this._PHI = PHI(this.x);
		var beta,res,PHIT;
		for (var i = 0; i < 10; i++){
			this._A = diag(this._alpha);
			this._SIGMA = getInverseMatrix(addMatrix(matrixMultiplication(1/this._sigma/this._sigma,matrixMultiplication(matrixT(this._PHI), this._PHI)), this._A));
			var_dump(this._SIGMA);
			beta = 1/this._sigma/this._sigma;
			var_dump('beta='+beta);
			PHIT = matrixT(this._PHI);
			var_dump("PHIT=",PHIT);
			res  = matrixMultiplication(matrixMultiplication(this._SIGMA,PHIT),matrixT(this.t));
			var_dump('res='+res);
			this._u = matrixMultiplication(beta, res);
			var_dump(this._u);
			this._gama = Array();
			for (var i = 0; i < this._SIGMA.length; i++){
				this._gama.push(1 - this._alpha[i] * this._SIGMA[i][i]);
			}
			var_dump(this._gama);
			for (var i = 0; i < this._alpha.length; i++){
				this._alpha[i] = this._gama[i] / this._u[i] / this._u[i];
			}
			this._sigma = F_Norm2(decMatrix(matrixT(this.t), matrixMultiplication(this._PHI, this._u))) / (this.t.length - E(this._gama, this.t.length));
			var_dump(this._alpha);
		}
	};
	this.forecast = function(x){};
}






















