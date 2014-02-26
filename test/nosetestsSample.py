def test_b():
    assert 'b' == 'b'
    print 'So what grade we got last year?'

class TestExampleTwo:
    def test_c(self):
        assert 'c' == 'c'

def multiply(a, b):
  """
  'multiply' multiplies two numbers and returns the result.

  >>> multiply(5, 10)
  50
  >>> multiply(-1, 1)
  -1
  >>> multiply(0.5, 1.5)
  0.75
  """
  return a*b

if __name__=='__main__':
    test_b()
