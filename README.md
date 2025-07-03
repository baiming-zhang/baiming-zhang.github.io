<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Baiming Zhang's GitHub Desktop</title>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 0;
            color: #333;
            background-color: #f9f9f9;
            line-height: 1.6;
        }
        .header {
            background-image: url('background.jpg');
            background-size: cover;
            background-position: center;
            color: white;
            padding: 30px 0;
            text-align: center;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }
        .container {
            max-width: 1000px;
            margin: 0 auto;
            padding: 20px;
        }
        nav {
            background-color: #0d47a1;
            padding: 15px 0;
        }
        nav a {
            color: white;
            text-decoration: none;
            margin: 0 20px;
            font-weight: 500;
            font-size: 16px;
            transition: color 0.3s;
        }
        nav a:hover {
            color: #bbdefb;
        }
        .profile {
            display: flex;
            flex-direction: column;
            align-items: center;
            text-align: center;
            padding: 40px 0;
            background-color: #fff;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
            margin-bottom: 30px;
        }
        .profile-photo {
            width: 147px;
            height: 206px;
            object-fit: cover;
            margin-bottom: 20px;
        }
        h1 {
            font-size: 2.5em;
            margin-bottom: 5px;
            color: #0d47a1;
        }
        h2 {
            font-size: 1.8em;
            color: #1a73e8;
            margin-top: 30px;
            border-bottom: 2px solid #e0e0e0;
            padding-bottom: 10px;
        }
        h3 {
            font-size: 1.4em;
            color: #2162b8; /* 天蓝色 */
        }
        .section {
            margin-bottom: 40px;
        }
        ul {
            padding-left: 20px;
        }
        li {
            margin-bottom: 15px;
        }
        .highlight {
            font-weight: bold;
            color: #d32f2f;
        }
        .timeline-item {
            margin-bottom: 25px;
        }
        footer {
            background-color: #333;
            color: white;
            text-align: center;
            padding: 20px 0;
            margin-top: 50px;
            box-shadow: 0 -2px 5px rgba(0, 0, 0, 0.2); /* 添加底部阴影 */
            width: 100%; /* 确保footer宽度为100% */
        }
        footer a {
            color: #bbdefb;
            text-decoration: none;
        }
        footer a:hover {
            text-decoration: underline;
        }
        .contact-btn {
            display: inline-block;
            background-color: #1a73e8;
            color: white;
            padding: 10px 20px;
            border-radius: 4px;
            text-decoration: none;
            margin-top: 15px;
            transition: background-color 0.3s;
        }
        .contact-btn:hover {
            background-color: #0d47a1;
        }
        @media (max-width: 768px) {
            .container {
                padding: 15px;
            }
            nav a {
                margin: 0 10px;
                font-size: 14px;
            }
            h1 {
                font-size: 2em;
            }
            h2 {
                font-size: 1.5em;
            }
        }
    </style>
</head>
<body>
    <header>
        <div class="container">
            <h1>🚀 Welcome to Baiming Zhang's Homepage!</h1>
        </div>
    </header>
    <nav>
        <div class="container">
            <a href="#home">Home</a>
            <a href="#Publications">Publications</a>
            <a href="#academic">Academic Experience</a>
            <a href="#awards">Awards & Honors</a>
            <a href="#boasts">Boasts</a>
            <a href="#contact">Contact</a>
        </div>
    </nav>
    <div class="container">
        <div class="profile" id="home">
            <img src="photo.jpg" alt="Baiming Zhang" class="profile-photo">
            <h2>Baiming Zhang</h2>
            <p>I'm a <span class="highlight">🎓 junior undergraduate</span> majoring in Engineering Mechanics at School of Aeronautics and Astronautics. 
            <p>I'm studying in Zhejiang University (ZJU) <span class="highlight">🚀 (Chu Kochen Honors College)</span> with a GPA of <span class="highlight">3.93/4.0</span> up to June 29, 2025.</p>
            <p>I'm currently conducting research in <span class="highlight">🔬 Quantum Computational Fluid Dynamics</span>.</p>
        </div>
        <div class="section" id="Publications">
            <h2>📰 Publications</h2>
            <p>At present, I have not yet published any papers. However, I am actively working on drafting the primary authorship of one paper, while also conceptualizing another that will also carry my primary authorship.</p>
            <ul>
                <li><strong>1 Drafting:</strong> Working on a paper as the primary author.</li>
                <li><strong>1 Conceptualizing:</strong> Developing ideas for another paper with primary authorship.</li>
            </ul>
            <p>I am excited about these projects and look forward to sharing more updates in the future.</p>
        </div>
        <div class="section" id="academic">
            <h2>📚 Academic and Research Experience</h2>
            <div class="timeline-item">
                <h3>1. ⚛️ Quantum States and Devices Research</h3>
                <ul>
                    <li>
                        <span class="highlight">🌟 October 2023 - June 2024:</span> <br>
                        Under the guidance of Professor Dawei Wang, Dr. Xingqi Xu and Dr. Jiefei Wang, I engaged in quantum optics research focusing on Rydberg atoms and topological insulators.
                    </li>
                    <li>
                        <span class="highlight">🧲 June 2024 - December 2024:</span> <br>
                        I explored magnetic field precision measurement and transitioned into studying Fock states and the 4D Quantum Hall Effect.
                    </li>
                </ul>
            </div>
            <div class="timeline-item">
                <h3>2. 🌊 Fluid Mechanics and Quantum Computing</h3>
                <ul>
                    <li>
                        <p><span class="highlight">   📚 January 2025 - April 2025:</span><br>
                        Under Professor Shiying Xiong's mentorship, I shifted my focus to fluid dynamics and quantum computing within my major. I also contributed to book writing related to my advisor's research.</p>
                    </li>
                </ul>
            </div>
            <div class="timeline-item">
                <h3>3. 🧠 Machine Learning and Quantum Hamilton Neural Networks (QHNN)</h3>
                <ul>
                    <li>
                        <span class="highlight">   🤖 April 2025 - June 2025:</span> <br>
                        I delved into neural networks, including Physics-Informed Neural Networks (PINN), Hamilton Neural Networks (HNN), and smaller networks like Taylor-Net.
                    </li>
                    <li>
                        <span class="highlight">   💡 June 2025 - Present:</span> <br>
                        Leveraging my quantum computing background and neural network expertise, I've been developing an original project called Quantum Hamiltonian Neural Networks (QHNN).
                    </li>
                </ul>
            </div>
        </div>
        <div class="section" id="awards">
            <h2>🏅 Awards and Honors</h2>
            <ol>
                <li>
                    <span class="highlight">🥇 International College Student Engineering Mechanics Competition</span> <br>
                    <strong>World Gold Medal (Top 0.1%)</strong>
                </li>
                <li>
                    <span class="highlight">🏆 Chinese National College Student Zhou Peiyuan Mechanics Competition</span> <br>
                    <strong>First Prize (Top 1%)</strong>
                </li>
                <li>
                    <span class="highlight">📜 Government Scholarship</span> <br>
                    <strong>(Top 10%)</strong>
                </li>
            </ol>
        </div>
        <div class="section" id="boasts">
            <h2>😊 Boasts</h2>
            <ol>
                <li>
                    <span class="highlight">🎉 Featured on Zhejiang University's School of Aeronautics and Astronautics official website</span> as the first student in our school's history to win the <strong>World Champion Gold Medal</strong> at the International College Student Engineering Mechanics Competition.
                </li>
                <li>
                    <span class="highlight">💡 Pursuing a breakthrough in cutting-edge quantum neural networks!</span> If you share this interest, I'd be delighted to connect with you. Feel free to reach out via email with any suggestions or inquiries.
                </li>
            </ol>
        </div>
    </div>
    <footer id="contact">
        <div class="container">
            <p>© 2025 Baiming Zhang. All rights reserved.</p>
            <p>Contact with me: <a href="mailto:baimingzhang@zju.edu.cn">baimingzhang@zju.edu.cn</a> | <a href="https://github.com/baiming-zhang">GitHub</a></p>
        </div>
    </footer>
</body>
</html>