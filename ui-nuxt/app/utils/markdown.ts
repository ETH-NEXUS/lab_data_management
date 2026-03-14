import DOMPurify from 'dompurify'
import { marked } from 'marked'

marked.setOptions({
  breaks: true,
  gfm: true,
})

const ALLOWED_TAGS = [
  'p',
  'br',
  'strong',
  'em',
  'u',
  's',
  'del',
  'ins',
  'h1',
  'h2',
  'h3',
  'h4',
  'h5',
  'h6',
  'ol',
  'ul',
  'li',
  'dl',
  'dt',
  'dd',
  'blockquote',
  'code',
  'pre',
  'span',
  'div',
  'a',
  'img',
  'table',
  'thead',
  'tbody',
  'tr',
  'th',
  'td',
  'hr',
  'sup',
  'sub',
]

export function renderMarkdown(content: string): string {
  if (!content) return ''

  try {
    const html = marked.parse(content, { async: false }) as string

    return DOMPurify.sanitize(html, {
      ALLOWED_TAGS,
      ALLOWED_ATTR: ['href', 'title', 'src', 'alt', 'class', 'align', 'width', 'height'],
    })
  } catch (e) {
    console.error('Markdown rendering failed:', e)
    return content
  }
}
