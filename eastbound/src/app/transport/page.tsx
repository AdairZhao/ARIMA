import type { Metadata } from "next";
import Link from "next/link";

export const metadata: Metadata = {
  title: "Transport",
  description:
    "High-speed rail, metro, DiDi, and airport transfers for travelers in China.",
};

const blocks = [
  {
    title: "High-speed rail",
    body: "Book early around national holidays. Trip.com is easier in English; 12306 is official. Arrive 45+ minutes early for security and real-name checks. Keep your passport out — it is your ticket ID.",
  },
  {
    title: "City metro",
    body: "Most hub cities have excellent metro networks. Add a transit QR code inside Alipay or the local metro mini-program. Watch the last-train times; DiDi fills late gaps.",
  },
  {
    title: "DiDi",
    body: "China’s main ride-hail app. English mode exists; destination pins are more reliable than typed addresses. Confirm the plate before you board.",
  },
  {
    title: "Airports",
    body: "Shanghai Maglev is optional, not mandatory. Beijing and Chengdu airports connect by metro or official taxi desks. Screenshot your hotel’s Chinese address for the driver.",
  },
];

export default function TransportPage() {
  return (
    <div className="wrap section--tight">
      <header className="page-hero">
        <p className="eyebrow">Transport</p>
        <h1 className="prose-title">How to move without Mandarin</h1>
        <p>
          Trains and metros are the backbone. Learn the passport check rhythm
          once — it repeats everywhere.
        </p>
      </header>

      <div className="check-list">
        {blocks.map((b) => (
          <article key={b.title} className="check-item">
            <h2 style={{ fontFamily: "var(--font-display)", fontSize: "1.45rem", margin: 0 }}>
              {b.title}
            </h2>
            <p className="muted" style={{ margin: 0 }}>
              {b.body}
            </p>
          </article>
        ))}
      </div>

      <p style={{ marginTop: "2rem" }}>
        <Link href="/essentials">Also read: payments &amp; apps →</Link>
      </p>
    </div>
  );
}
