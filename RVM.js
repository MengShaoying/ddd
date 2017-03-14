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
	var h = '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" height="'
		+ (k*dy+5) + '" width="' + k*dx +'">';
	for (var i = 0; i < x.length; i++){
		h += '<circle cx="'+k*(x[i]-min(x))+'" cy="'
			+ k*(dy-y[i]+min(y)) + '" r="2" fill="black"/>';
	}
	h += '</svg>';
	var p = document.createElement('div');
	p.innerHTML = h;
	p.setAttribute("class","graph");
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
		return Math.exp((0-r2)/(2*sigma*sigma)) /
			Math.sqrt(2*Math.PI*sigma);
	} else {
		for (var i = 0; i < x.length; i++)
			r2 += (x[i] - mu[i]) * (x[i] - mu[i]);
		return Math.exp((0-r2)/(2*sigma*sigma)) / 
			Math.pow(Math.sqrt(2*Math.PI*sigma),x.length);
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
			html += matrix[i][j].toFixed(2) + '&nbsp;';
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
	if ('number' == typeof(a)){
		for (var i = 0; i < b.length; i++)
			for (var j = 0; j < b[i].length; j++)
				b[i][j] *= a;
		return b;
	}
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

function array2matrix(t){
	var a0 = Array();
	for (var i = 0; i < t.length; i++){
		a0.push(t[i]);
	}
	return [a0];
}

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
			s.push(gaussianFunction(x[i], x[j], 4));
		p.push(s);
	}
	return p;
}

function RVM(){
	this.x = Array();
	this.t = Array();
	this._alpha = Array();/* 一维数组 */
	this._sigma = 0;      /* 实数 */
	this._PHI   = Array();/* 矩阵 */
	this._PHIT  = Array();/* 矩阵 */
	this._A     = Array();/* 矩阵 */
	this._SIGMA = Array();/* 矩阵 */
	this._u     = Array();/* 矩阵 */
	this._gama  = Array();/* 一维数组 */
	this.setTrainXT = function(x,t){
		this.x  = x;
		this.t  = t;
	};
	this.initAlphaSigma = function(){
		this._alpha = Array();
		for (var i = 0; i < this.x.length + 1; i++)
			this._alpha.push(1 + Math.random());
		this._sigma = 1 + Math.random();
	};
	this.train = function(){
		this._PHI  = PHI(this.x);
		this._PHIT = matrixT(this._PHI);
		var t = array2matrix(this.t);
		var beta, ktxk, top, buttom;
		for (var i = 0; i < 100; i++){
			//var_dump('i=');
			//var_dump(i);
			this._A = diag(this._alpha);
			//var_dump('A=');
			//viewMatrix(this._A);
			/* 计算SIGMA */
			beta = 1 / this._sigma / this._sigma;
			//var_dump('beta=');
			//var_dump(beta);
			ktxk = matrixMultiplication(this._PHIT, this._PHI);
			ktxk = matrixMultiplication(beta, ktxk);
			ktxk = addMatrix(ktxk, this._A);
			this._SIGMA = getInverseMatrix(ktxk);
			//var_dump('SIGMA=');
			//viewMatrix(this._SIGMA);
			/* 计算this._u和this._gama */
			this._u = matrixMultiplication(beta, this._SIGMA);
			this._u = matrixMultiplication(this._u, this._PHIT);
			this._u = matrixMultiplication(this._u, matrixT(t));
			//var_dump('u=');
			//viewMatrix(this._u);
			/* 计算gama值 */
			this._gama = Array();
			for (var k = 0; k < this._alpha.length; k++){
				this._gama.push(
					1 - this._alpha[k] * this._SIGMA[k][k]
				);
			}
			//var_dump('gama=');
			//var_dump(this._gama);
			/* 计算新的alpha值 */
			for (var k = 0; k < this._alpha.length; k++){
				this._alpha[k] = this._gama[k] /
					(this._u[k][0] * this._u[k][0]);
			}
			//var_dump('new alpha=');
			//var_dump(this._alpha);
			/* 计算新的sigma值 */
			top = matrixMultiplication(this._PHI, this._u);
			top = matrixMultiplication(-1, top);
			top = addMatrix(matrixT(t), top);
			top = F_Norm2(top);
			buttom = 0;
			for (var k = 0; k < this._gama.length; k++){
				buttom += this._gama[k];
			}
			//var_dump('N=');
			//var_dump(this.x.length);
			//var_dump('top=');
			//var_dump(top);
			//var_dump('buttom=');
			//var_dump(buttom);
			buttom = this.x.length - buttom;
			if (buttom < 0){
				//var_dump('i=');
				//var_dump(i);
				//var_dump('N-gama[i]*SIGMA[ii]<0, end.');
				break;
			}
			buttom = buttom < 0 ? -buttom : buttom;
			this._sigma = Math.sqrt(top / buttom);
			//var_dump('new sigma=');
			//var_dump(this._sigma);
			//var_dump('*******************************************');
		}
	};
	this.forecast = function(x){
		/* 计算phi */
		var k = [1];
		for (var i = 0; i < this.x.length; i++){
			k.push(gaussianFunction(x, this.x[i], 4));
		}
		//var_dump(k);
		k = array2matrix(k);
		k = matrixT(k);
		//viewMatrix(k);
		/* 计算预测结果:uT*phi(x) */
		var r = matrixMultiplication(matrixT(this._u),k);
		return r[0][0];
	};
}






















