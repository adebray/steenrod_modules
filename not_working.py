#!/usr/bin/env python3

# code for computing resolutions over the truncation A(1) of the Steenrod algebra at p = 2
# 21 May 2017

import copy
import random


# elements of an A(1)-module
class ModElem:
	sq1: str = None # tag of Sq1 of this element. None => zero
	sq2: str = None # tag of Sq2 of this element. None => zero
	degree: int = None # degree of the element
	in_kernel: bool = None # if we're calculating a resolution, is this in the kernel?
	is_hit_by: str = None # if we're calculating a resolution, what's the tag of the thing that hits this?
	duplicate: str = None # if we're calculating a resolution, does this map to the same place as anything?

	def __init__(self, degree, sq1, sq2):
		self.sq1 = sq1
		self.sq2 = sq2
		self.degree = degree
		self.in_kernel = None
		self.is_hit_by = None
		self.duplicate = None

	def __str__(self):
		return ('ModElem %d -> %s %s; ker: %r hit: %s dup: %s' %
			(self.degree, self.sq1, self.sq2, self.in_kernel, self.is_hit_by, self.duplicate))
	def __repr__(self):
		return str(self)

# the zero element
O = ModElem(0, None, None)
O.sq1 = O
O.sq2 = O

# not quite a defaultdict: doesn't add the missing entry
# https://stackoverflow.com/a/17956989
class default_to_zero(dict):
	def __missing__(self, key):
		return O
	
# An A(1)-module is a list of ModElems (I'm not doing any validation, oops)
# Here are some examples
A1 = {
	'A0': ModElem(0, 'A1', 'A2'),
	'A1': ModElem(1, None, 'B3'),
	'A2': ModElem(2, 'A3', 'B4'),
	'A3': ModElem(3, None, 'B5'),
	'B3': ModElem(3, 'B4', None),
	'B4': ModElem(4, None, 'B6'),
	'B5': ModElem(5, 'B6', None),
	'B6': ModElem(6, None, None)
}

A0 = {
	'A0': ModElem(0, 'A1', None),
	'A1': ModElem(1, None, None)
}

C2 = {
	'C0': ModElem(0, None, 'C2'),
	'C2': ModElem(0, None, None)
}

F2 = {
	'F': ModElem(0, None, None)
}

# we disambiguate different copies of A1 using 'tags,' which are strings that are long enough to probably be
# unique in each run
def new_tag() -> str:
	return '-' + hex(random.randrange(1e8))[2:]

# modifies x so that sq^1x and sq^2x have tags appended to them
def append_tag(x: ModElem, tag: str):
	if x.sq1 is not None: x.sq1 += tag
	if x.sq2 is not None: x.sq2 += tag
	return x

# produces a deep copy of A1 whose labels have been tagged. Returns the tag as the second argument
def new_A1_copy() -> (list, str):
	tag = new_tag()
	A1_copy = copy.deepcopy(A1)
	to_return = dict((name + tag, append_tag(val, tag)) for (name, val) in A1_copy.items())
	return to_return, tag

# suspend the A(1)-module M n times.
# the original still exists if duplicate=True
def suspend(M: dict, n=1, duplicate=True) -> dict:
	to_return = copy.deepcopy(M) if duplicate else M
	for tag in to_return:
		to_return[tag].degree += n
	return to_return

# args should be A(1)-modules, i.e. dicts of ModElems
def direct_sum(*args):
	to_return = default_to_zero()
	for M in args:
		to_return.update(M)
	return to_return

# returns an element of M of minimal degree that's not already in the image of something
def smallest_yet_to_hit(M):
	argmin_so_far = None
	argmin_tag = None
	val = 1e10
	for tag, m in M.items():
		#print(tag)
		if (m.in_kernel or m.duplicate) and not m.is_hit_by:
			next_val = m.degree
			if next_val < val:
				argmin_so_far = m
				argmin_tag = tag
				val = next_val
	return argmin_tag, argmin_so_far

# maps sq^i to the elements of s. Removes all zeros
def apply_sq(i: int, s: set):
	sq = lambda x: x.sq1 if i = 1 else lambda x: x.sq2
	return set(sq(x) for x in s if x != O)

# checks whether 'to' is zero, and updates 'from' accordingly
# TODO: eventually, replace in_kernel with image, which could be 0
# returns to.duplicate, so that later calls to elem_update can use it
# note: image is a set of what fr_tag hits and its duplicates
def elem_update(M, M_next, image, fr_tag)
	to_return = image
	fr = M_next[fr_tag]
	if not dupe_check:
		fr.in_kernel = (to == O)
	for to in image:
		if to != O:
			print(fr_tag, '--> ', to)
			if to.is_hit_by:
				fr.duplicate = to.is_hit_by
				M_next[to.is_hit_by].duplicate = fr_tag
				print('\033[33mfound duplicate: %s _hit by_ %s _and_ %s\033[0m' % (to, to.is_hit_by, fr_tag))
			to.is_hit_by = fr_tag
			if to.duplicate:
				M[to.duplicate].is_hit_by = fr_tag
				to_return.add(to.duplicate)
				print('\033[33mhit duplicate: %s _and_ %s _hit by_ %s\033[0m' % (to, to.duplicate, fr_tag))
	return to_return
		
# TODO
# current issue: if x maps to two things, sq^1x and sq^2x must map to the images of those things
# potential way to solve this: elem_update can return an 'extra' value that we also have to hit

# how about: elem_update takes a list of (one or two) targets?
# in kernel if all are O.
# I should be able to prevent the number of duplicates from spiraling out of control

# computes the image given by adding an A(1) to M_next
def update_info(M, m_tag, M_next, new_A1_tag):
	m0 = M[m_tag]

	dupes_A0 = elem_update(M, M_next, set(m0), 'A0' + new_A1_tag)


	elem_update(M, M_next, M[m0.sq1], 'A1' + new_A1_tag)
	elem_update(M, M_next, M[m0.sq2], 'A2' + new_A1_tag)
	elem_update(M, M_next, M[M[m0.sq2].sq1], 'A3' + new_A1_tag)
	elem_update(M, M_next, M[M[m0.sq1].sq2], 'B3' + new_A1_tag)
	elem_update(M, M_next, M[M[M[m0.sq1].sq2].sq1], 'B4' + new_A1_tag)
	elem_update(M, M_next, M[M[M[m0.sq2].sq1].sq2], 'B5' + new_A1_tag)
	elem_update(M, M_next, M[M[M[M[m0.sq1].sq2].sq1].sq2], 'B6' + new_A1_tag)

# computes the next step in a minimal projective resolution, where M is the previous step
# returns an A(1)-module
def minimal_resolution_step(M):
	M_next = default_to_zero()
	print('next step:')
	while True:
		# 1. find a lowest-degree element of M that's not already hit or 
		m_tag, m = smallest_yet_to_hit(M)
		if m_tag is None: return M_next
		#print(m, m_tag)

		# 2. hit it with an A(1)
		new_A1, new_A1_tag = new_A1_copy()
		new_A1 = suspend(new_A1, m.degree, duplicate=False)
		print('\tΣ^%d 𝓐(1)' % m.degree)
	
		# 3. compute the image
		M_next = direct_sum(M_next, new_A1)
		update_info(M, m_tag, M_next, new_A1_tag)

def minimal_resolution(M):
	for key in M:
		M[key].in_kernel = True
	to_return = list(M)
	while True:
		try:
			M_next = minimal_resolution_step(M)
			#for k, v in M_next.items():
			#	print(k, v)
			#print(list((k, str(v)) for k, v in M_next.items()))
			input()
			to_return.append(M_next)
			M = M_next
		except KeyboardInterrupt:
			#print(to_return)
			break

def main():
	M = default_to_zero(A0)
	minimal_resolution(M)

if __name__ == '__main__':
	main()
