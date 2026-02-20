"""
Biotechnology Knowledge Base - Flask Web Application
"""

from flask import Flask, render_template, request, jsonify
import sqlite3
import os

app = Flask(__name__)

DB_PATH = "biotech_kb.db"

def get_db_connection():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn

@app.route('/')
def index():
    """Home page"""
    conn = get_db_connection()
    cursor = conn.cursor()
    
    # Get total terms
    cursor.execute('SELECT COUNT(*) as count FROM terms')
    total_terms = cursor.fetchone()[0]
    
    # Get categories with counts
    cursor.execute('SELECT category, COUNT(*) as count FROM terms GROUP BY category ORDER BY count DESC')
    categories = [dict(row) for row in cursor.fetchall()]
    
    conn.close()
    
    return render_template('kb_index.html', 
                         total_terms=total_terms, 
                         categories=categories)

@app.route('/search')
def search():
    """API endpoint for searching terms"""
    query = request.args.get('q', '').strip()
    category = request.args.get('category', '')
    limit = int(request.args.get('limit', 20))
    
    conn = get_db_connection()
    cursor = conn.cursor()
    
    if query and category:
        cursor.execute('''
            SELECT * FROM terms 
            WHERE (term LIKE ? OR definition LIKE ?) AND category = ?
            LIMIT ?
        ''', (f'%{query}%', f'%{query}%', category, limit))
    elif query:
        cursor.execute('''
            SELECT * FROM terms 
            WHERE term LIKE ? OR definition LIKE ?
            LIMIT ?
        ''', (f'%{query}%', f'%{query}%', limit))
    elif category:
        cursor.execute('SELECT * FROM terms WHERE category = ? LIMIT ?', (category, limit))
    else:
        cursor.execute('SELECT * FROM terms LIMIT ?', (limit,))
    
    results = [dict(row) for row in cursor.fetchall()]
    conn.close()
    
    return jsonify(results)

@app.route('/term/<int:term_id>')
def get_term(term_id):
    """Get single term details"""
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM terms WHERE id = ?', (term_id,))
    term = dict(cursor.fetchone()) if cursor.fetchone() else None
    conn.close()
    
    if term:
        return jsonify(term)
    return jsonify({'error': 'Term not found'}), 404

@app.route('/categories')
def get_categories():
    """Get all categories"""
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute('SELECT category, COUNT(*) as count FROM terms GROUP BY category ORDER BY count DESC')
    categories = [dict(row) for row in cursor.fetchall()]
    conn.close()
    
    return jsonify(categories)

@app.route('/stats')
def stats():
    """Get database statistics"""
    conn = get_db_connection()
    cursor = conn.cursor()
    
    cursor.execute('SELECT COUNT(*) as total FROM terms')
    total = cursor.fetchone()[0]
    
    cursor.execute('SELECT category, COUNT(*) as count FROM terms GROUP BY category')
    by_category = [dict(row) for row in cursor.fetchall()]
    
    cursor.execute('SELECT term_type, COUNT(*) as count FROM terms GROUP BY term_type')
    by_type = [dict(row) for row in cursor.fetchall()]
    
    conn.close()
    
    return jsonify({
        'total_terms': total,
        'by_category': by_category,
        'by_type': by_type
    })

if __name__ == '__main__':
    print("Starting Biotechnology Knowledge Base...")
    print(f"Database: {DB_PATH}")
    app.run(debug=True, host='0.0.0.0', port=5001)
