import type { Metadata } from "next";
import Link from "next/link";
import { planChecklist } from "@/content/checklist";

export const metadata: Metadata = {
  title: "Plan checklist",
  description:
    "Visa, payments, eSIM, apps, and train booking — the pre-trip checklist for visiting China.",
};

export default function PlanPage() {
  return (
    <div className="wrap section--tight">
      <header className="page-hero">
        <p className="eyebrow">Plan</p>
        <h1 className="prose-title">Clear the gates before you fly</h1>
        <p>
          Most first trips fail on logistics, not on sights. Work this list in
          order — then pick cities with confidence.
        </p>
      </header>

      <div className="check-list">
        {planChecklist.map((item) => (
          <article key={item.id} className="check-item">
            <div className="check-item__top">
              <h2 style={{ fontFamily: "var(--font-display)", fontSize: "1.45rem", margin: 0 }}>
                {item.title}
              </h2>
              <span className="timing">{item.timing}</span>
            </div>
            <p className="muted" style={{ margin: 0 }}>
              {item.detail}
            </p>
          </article>
        ))}
      </div>

      <div className="cta-band" style={{ marginTop: "3rem" }}>
        <h2>Next: choose your hubs</h2>
        <p>Four starter cities with see / eat / move notes in English.</p>
        <Link href="/destinations" className="btn btn--primary">
          Browse destinations
        </Link>
      </div>
    </div>
  );
}
